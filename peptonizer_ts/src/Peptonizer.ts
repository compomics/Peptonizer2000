import {GridSearchProgressListener} from "./GridSearchProgressListener.ts";
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
     * @return Mapping between NCBI taxon IDs (integer, > 0) and probabilities (float in [0, 1]).
     */
    async peptonize(
        peptidesScores: Map<string, number>,
        peptidesCounts: Map<string, number>,
        alphas: number[],
        betas: number[],
        priors: number[],
        progressListener: GridSearchProgressListener
    ): Promise<PeptonizerResult[]> {
        const workerPool = new WorkerPool();

        const peptonizerResults: PeptonizerResult[] = [];

        try {
            const taxonWeightsCsv = await workerPool.performTaxaWeighing(peptidesScores, peptidesCounts);
            const generatedGraph = await workerPool.generateGraph(taxonWeightsCsv);

            for (const alpha of alphas) {
                for (const beta of betas) {
                    for (const prior of priors) {
                        const result = await workerPool.executePepgm(
                            generatedGraph,
                            alpha,
                            beta,
                            prior,
                            progressListener
                        );

                        peptonizerResults.push(result);
                    }
                }
            }
        } catch (err) {
            console.error(err);
        }

        return peptonizerResults
    }
}

export { Peptonizer };
export type { PeptonizerResult }