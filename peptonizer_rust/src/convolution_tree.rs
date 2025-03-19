use rustfft::{FftPlanner, num_complex::Complex, num_traits::Zero};
use crate::array_utils::{normalize, log_normalize};

#[derive(Debug, Clone)]
struct CTNode {
    joint_above: Vec<f64>,
    left_parent: Option<Box<CTNode>>,
    right_parent: Option<Box<CTNode>>,
    likelihood_below: Option<Vec<f64>>,
}

impl CTNode {
    /// Create a new CTNode with a given probability vector.
    fn new(mut joint_above: Vec<f64>) -> Self {
        normalize(&mut joint_above);
        Self {
            joint_above: joint_above,
            left_parent: None,
            right_parent: None,
            likelihood_below: None,
        }
    }

    /// Creates a count node by convolving two parent nodes.
    fn create_count_node(lhs: CTNode, rhs: CTNode) -> CTNode {
        let joint_above = fft_convolve(&lhs.joint_above, &rhs.joint_above);
        let mut node = CTNode::new(joint_above);
        node.left_parent = Some(Box::new(lhs));
        node.right_parent = Some(Box::new(rhs));
        node
    }

    /// Computes the upward message by convolving with the given probability vector.
    fn message_up(&self, answer_size: usize, other_joint_vector: &Vec<f64>) -> Vec<f64> {
        let likelihood = self.likelihood_below.as_ref().expect("Likelihood below is None!");
        let starting_point = other_joint_vector.len() - 1;
        let result = fft_convolve(
            &other_joint_vector.iter().rev().cloned().collect::<Vec<f64>>(),
            likelihood,
        );
        let mut result = result[starting_point..starting_point + answer_size].to_vec();
        normalize(&mut result);
        result
    }

    /// Computes message passing upwards to the left parent.
    fn message_up_left(&self) -> Vec<f64> {
        let right_parent = self.right_parent.as_ref().expect("Right parent is missing!");
        let left_parent = self.left_parent.as_ref().expect("Left parent is missing!");
        self.message_up(left_parent.joint_above.len(), &right_parent.joint_above)
    }

    /// Computes message passing upwards to the right parent.
    fn message_up_right(&self) -> Vec<f64> {
        let left_parent = self.left_parent.as_ref().expect("Left parent is missing!");
        let right_parent = self.right_parent.as_ref().expect("Right parent is missing!");
        self.message_up(right_parent.joint_above.len(), &left_parent.joint_above)
    }

    /// Computes posterior probability.
    fn posterior(&self) -> Vec<f64> {
        let likelihood = self.likelihood_below.as_ref().expect("Likelihood below is None!");
        let mut result = self.joint_above.iter().zip(likelihood.iter()).map(|(a, b)| a * b).collect();
        normalize(&mut result);
        result
    }

    /// Returns messages up.
    fn messages_up(&self) -> Vec<f64> {
        self.likelihood_below.clone().expect("Likelihood below is None!")
    }
}

#[derive(Debug)]
pub struct ConvolutionTree {
    n_to_shared_likelihoods: Vec<f64>,
    log_length: usize,
    all_layers: Vec<Vec<CTNode>>,
    last_node: Option<CTNode>,
    protein_layer: Vec<CTNode>,
    n_proteins: usize
}

impl ConvolutionTree {
    /// Constructs a ConvolutionTree.
    pub fn new(n_to_shared_likelihoods: Vec<f64>, proteins: Vec<Vec<f64>>) -> Self {
        let log_length = (proteins.len() as f64).log2().ceil() as usize;
        let mut tree = ConvolutionTree {
            n_to_shared_likelihoods,
            log_length,
            all_layers: Vec::new(),
            last_node: None,
            protein_layer: Vec::new(),
            n_proteins: proteins.len()
        };

        tree.build_first_layer(proteins);
        tree.build_remaining_layers();
        tree.propagate_backward();

        tree
    }

    /// Builds the first layer of the tree (protein nodes).
    fn build_first_layer(&mut self, proteins: Vec<Vec<f64>>) {
        let mut layer = proteins.into_iter().map(CTNode::new).collect::<Vec<CTNode>>();

        // Pad with dummy nodes to make the length a power of 2
        while layer.len() < 2usize.pow(self.log_length as u32) {
            layer.push(CTNode::new(vec![1.0, 0.0]));
        }

        self.all_layers.push(layer);
    }

    /// Builds the remaining layers using count nodes.
    fn build_remaining_layers(&mut self) {
        for _ in 0..self.log_length {
            let most_recent_layer = self.all_layers.last().unwrap();
            let mut new_layer = Vec::new();

            for i in (0..most_recent_layer.len()).step_by(2) {
                let left = most_recent_layer[i].clone();
                let right = most_recent_layer[i + 1].clone();
                new_layer.push(CTNode::create_count_node(left, right));
            }

            self.all_layers.push(new_layer);
        }

        let mut final_node = self.all_layers.last().unwrap()[0].clone();
        let mut likelihood_below = self.n_to_shared_likelihoods.clone();
        normalize(&mut likelihood_below);
        final_node.likelihood_below = Some(likelihood_below);

        self.all_layers.last_mut().unwrap()[0] = final_node.clone();
        self.last_node = Some(final_node);
    }

    /// Propagates messages backwards through the tree.
    fn propagate_backward(&mut self) {
        for L in (1..=self.log_length).rev() {
            let layer = &mut self.all_layers[L];

            for node in layer.iter_mut() {
                let likelihood_below_left = Some(node.message_up_left());
                let likelihood_below_right = Some(node.message_up_right());

                node.left_parent.as_mut().unwrap().likelihood_below = likelihood_below_left;
                node.right_parent.as_mut().unwrap().likelihood_below = likelihood_below_right;
            }
        }

        self.protein_layer = self.all_layers[0].clone();
    }

    fn posterior_for_variable(&self, prot_idx: usize) -> Vec<f64> {
        self.protein_layer[prot_idx].posterior()
    }

    pub fn message_to_variable(&self, prot_idx: usize) -> Vec<f64> {
        self.protein_layer[prot_idx].messages_up()
    }

    pub fn message_to_shared_likelihood(&self) -> Vec<f64> {
        let last_node = self.last_node.as_ref().expect("Last node is missing");

        // Extract the required range
        last_node.joint_above[..=self.n_proteins].to_vec()
    }
}

/// Performs FFT-based convolution using `rustfft`.
fn fft_convolve(a: &Vec<f64>, b: &Vec<f64>) -> Vec<f64> {
    let len = a.len() + b.len() - 1;
    let fft_size = len.next_power_of_two();

    let mut planner = FftPlanner::new();
    let fft = planner.plan_fft_forward(fft_size);
    let ifft = planner.plan_fft_inverse(fft_size);

    let mut a_complex: Vec<Complex<f64>> = a.iter().map(|&x| Complex::new(x, 0.0)).collect();
    let mut b_complex: Vec<Complex<f64>> = b.iter().map(|&x| Complex::new(x, 0.0)).collect();
    a_complex.resize(fft_size, Complex::zero());
    b_complex.resize(fft_size, Complex::zero());

    fft.process(&mut a_complex);
    fft.process(&mut b_complex);

    let mut result_complex: Vec<Complex<f64>> = a_complex.iter().zip(&b_complex).map(|(x, y)| x * y).collect();
    ifft.process(&mut result_complex);

    result_complex.iter().take(len).map(|c| c.re / fft_size as f64).collect()
}