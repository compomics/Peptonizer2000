import async from "async";
import { GridSearchProgressListener } from "./GridSearchProgressListener.ts";
import ExecutePepGMWorker from './workers/ExecutePepGMWorker.ts?worker&inline';

type BeliefPropagationParameters = {
    alpha: number,
    beta: number,
    prior: number
}

type BeliefPropagationTask = {
    graph: string,
    parameters: BeliefPropagationParameters
    progressListener: GridSearchProgressListener
}

// Map of organism name to confidence score computed by PepGM
type BeliefPropagationResult = Map<string, number>;

class GridSearchWorkerPool {
    static callbacks = new Map();
    static currentId = 0;

    // Maximum amount of workers that are allowed to be run in parallel (using multiple CPU cores).
    static numberOfWorkers = 1;
    static workers: Worker[] = [];

    static queue;

    // Properly initialize the array of workers and their logic required to handle their asynchronous nature.
    static {
        for (let i = 0; i < this.numberOfWorkers; i++) {
            const worker = new ExecutePepGMWorker();

            // Define what to do if a pepGM execution is done.
            worker.onmessage = (event) => {
                const {id, type, ...data} = JSON.parse(event.data);
                const [onSuccess, onError, progressListener] = this.callbacks.get(id);

                if (type === "result") {
                    this.callbacks.delete(id);
                    onSuccess(JSON.parse(data.result));
                } else if (type === "progress") {
                    if (data["progress_type"] === "graph") {
                        progressListener.graphsUpdated(
                            data["current_graph"],
                            data["total_graphs"],
                            id % this.numberOfWorkers
                        );
                    } else if (data["progress_type"] === "max_residual") {
                        progressListener.maxResidualUpdated(
                            data["max_residual"],
                            data["tolerance"],
                            id % this.numberOfWorkers
                        );
                    } else if (data["progress_type"] === "iterations") {
                        progressListener.iterationsUpdated(
                            data["current_iterations"],
                            data["total_iterations"],
                            id % this.numberOfWorkers
                        )
                    }
                } else if (type === "error") {
                    this.callbacks.delete(id);
                    onError(data.error);
                }
            };

            this.workers.push(worker);
        }

        this.queue = async.queue(
            (
                task: BeliefPropagationTask,
                callback: (x: BeliefPropagationResult | { error: Error }) => void
            ) => {
                (async () => {
                    // Retrieve worker from the pool.
                    const worker = this.workers.pop();

                    if (worker) {
                        try {
                            const result: BeliefPropagationResult = await new Promise((onSuccess, onError) => {
                                this.currentId = (this.currentId + 1) % Number.MAX_SAFE_INTEGER;
                                this.callbacks.set(
                                    this.currentId,
                                    [onSuccess, onError, task.progressListener]
                                );

                                const params = task.parameters;

                                // Notify the progress listener that we've started training a new graph with different
                                // parameters.
                                task.progressListener.gridUpdated(
                                    params,
                                    this.currentId % this.numberOfWorkers
                                );

                                // The progressListener itself cannot be sent to the worker!
                                worker.postMessage({
                                    graph: task.graph,
                                    alpha: task.parameters.alpha,
                                    beta: task.parameters.beta,
                                    prior: task.parameters.prior,
                                    id: this.currentId
                                });
                            });

                            // Add worker back to the pool.
                            this.workers.push(worker);

                            // We're done, return the results to the callback
                            callback(result);
                        } catch (err) {
                            callback({ error: err as Error })
                        }
                    } else {
                        callback({ error: new Error("No workers available in the queue!") })
                    }
                })();
            }, this.numberOfWorkers
        );
    }

    /**
     * Run the belief propagation algorithm on the provided graph for all possible combinations of the parameter
     * values.
     *
     * @param graph
     * @param alphas
     * @param betas
     * @param priors
     * @param progressListener
     * @returns {Promise}
     */
    static performGridSearch(
        graph: string,
        alphas: number[],
        betas: number[],
        priors: number[],
        progressListener: GridSearchProgressListener
    ): Promise<BeliefPropagationResult[]> {
        return new Promise((resolve, reject) => {
            const results: BeliefPropagationResult[] = [];

            // Loop through all parameter combinations.
            for (const alpha of alphas) {
                for (const beta of betas) {
                    for (const prior of priors) {
                        // Push the tasks to the queue.
                        this.queue.push(
                            {
                                graph,
                                parameters: {
                                    alpha,
                                    beta,
                                    prior
                                },
                                progressListener
                            },
                            // @ts-ignore
                            (result: BeliefPropagationResult | { error: Error }) => {
                                // Check if an error occurred. If so, stop running the grid search entirely.
                                // @ts-ignore
                                if (result["error"]) {
                                    // @ts-ignore
                                    reject(result.error);
                                    return;
                                }

                                // No error, continue and add the results to the output
                                results.push(result as BeliefPropagationResult)
                            }
                        );
                    }
                }
            }

            // send pepGM onSucces if all tasks are processed
            this.queue.drain(() => {
                resolve(results);
            });
        });
    }
}

export { GridSearchWorkerPool };
export type { BeliefPropagationResult, BeliefPropagationParameters };
