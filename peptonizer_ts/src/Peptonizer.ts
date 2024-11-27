import {PeptonizerParameterSet, PeptonizerProgressListener} from "./PeptonizerProgressListener.ts";
import {WorkerPool} from "./workers/WorkerPool.ts";

type PeptonizerResult = Map<string, number>;

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
     * @param workers The amount of Web Workers that can be spawned and used simultaneously to run the Peptonizer.
     * @return Mapping between NCBI taxon IDs (integer, > 0) and probabilities (float in [0, 1]).
     */
    async peptonize(
        peptidesScores: Map<string, number>,
        peptidesCounts: Map<string, number>,
        alphas: number[],
        betas: number[],
        priors: number[],
        progressListener: PeptonizerProgressListener,
        workers: number = 2
    ): Promise<PeptonizerResult> {
        const workerPool = new WorkerPool(workers);

        const [sequenceScoresCsv, taxonWeightsCsv] = await workerPool.performTaxaWeighing(peptidesScores, peptidesCounts);
        const generatedGraph = await workerPool.generateGraph(sequenceScoresCsv);

        // Compute the total amount of tasks and which tasks to perform.
        const parameterSets: PeptonizerParameterSet[] = [];

        for (const alpha of alphas) {
            for (const beta of betas) {
                for (const prior of priors) {
                    parameterSets.push({
                        alpha,
                        beta,
                        prior
                    });
                }
            }
        }

        // Notify any listeners that the Peptonizer did start running (and report which set of parameters will be tuned)
        progressListener?.peptonizerStarted(parameterSets.length, parameterSets);

        const pepgmPromises: Promise<PeptonizerResult>[] = [];

        for (const paramSet of parameterSets) {
            pepgmPromises.push(
                workerPool.executePepgm(generatedGraph, paramSet.alpha, paramSet.beta, paramSet.prior, progressListener)
            );
        }

        // Wait until all parameter sets have been tuned...
        const peptonizerResults = await Promise.all(pepgmPromises);

        // Now that we have all the results generated by the peptonizer, we need to figure out which one yields the best
        // results

        // First compute the clustered taxa weights
        const clusteredTaxaWeightsCsv = await workerPool.clusterTaxa(generatedGraph, taxonWeightsCsv);

        let bestGoodness = -1;
        let bestResult: PeptonizerResult | undefined;
        for (const result of peptonizerResults) {
            const goodness = await workerPool.computeGoodness(clusteredTaxaWeightsCsv, result);

            if (goodness > bestGoodness) {
                bestGoodness = goodness;
                bestResult = result;
            }
        }

        if (!bestResult) {
            throw new Error("No results found!");
        }

        progressListener?.peptonizerFinished();

        return bestResult;
    }
}

export { Peptonizer };
export type { PeptonizerResult }