import {GridSearchProgressListener} from "./GridSearchProgressListener.ts";
import {GraphGenerationWorkerPool} from "./GraphGenerationWorkerPool.ts";
import {BeliefPropagationResult, GridSearchWorkerPool} from "./GridSearchWorkerPool.ts";

class Peptonizer {
    /**
     * Start the peptonizer. This function takes in a PSM-file that has been read in earlier (so no file paths here). The
     * PSMS including their intensities are then used as input to the Peptonizer-pipeline. This pipeline will finally
     * return a Map in which NCBI taxon IDs are mapped onto their propabilities (as computed by the Peptonizer).
     *
     * @param peptidesScores Mapping between the peptide sequences that should be present in the peptonizer and a scoring
     * metric (derived from a prior search engine step) for each peptide.
     * @param peptidesCounts Mapping between peptide sequences and the amount of times they occur in the input.
     * @param alphas An array of possible values for the alpha parameter. This parameter indicates the probability that an
     * observed taxon also indicates the presence of a peptide.
     * @param betas An array of possible values for the beta parameter. This parameter indicates the probability of
     * detecting a peptide at random.
     * @param priors An array of possible values for the gamma (or prior) parameter. Gamma indicates the prior probability
     * of a taxon being present.
     * @param progressListener Is called everytime the progress of the belief propagation algorithm has been updated.
     * @return Mapping between NCBI taxon IDs (integer, > 0) and probabilities (float in [0, 1]).
     */
    async peptonize(
        peptidesScores: Map<string, number>,
        peptidesCounts: Map<string, number>,
        alphas: number[],
        betas: number[],
        priors: number[],
        progressListener: GridSearchProgressListener
    ): Promise<BeliefPropagationResult[]> {
        // First, we've got to generate the PepGM factor graph itself.
        const generatedGraph = await GraphGenerationWorkerPool.generatePepGmGraph(peptidesScores, peptidesCounts);

        // // Check if generating the graph worked properly
        // if (generatedGraph.error) {
        //     return {
        //         error: generatedGraph.error
        //     }
        // }

        // If we reach this point, the graphs have been generated properly, and we can start the belief propagation
        // algorithm for all possible parameter combinations.
        return GridSearchWorkerPool.performGridSearch(
            generatedGraph,
            alphas,
            betas,
            priors,
            progressListener
        );
    }
}

export { Peptonizer };