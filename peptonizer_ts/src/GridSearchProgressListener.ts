import { BeliefPropagationParameters } from "./GridSearchWorkerPool.ts";

interface GridSearchProgressListener {
    gridUpdated(params: BeliefPropagationParameters, workerId: number): void;
    graphsUpdated(currentGraph: number, totalGraphs: number, workerId: number): void;
    maxResidualUpdated(maxResidual: number, tolerance: number, workerId: number): void;
    iterationsUpdated(currentIteration: number, totalIterations: number, workerId: number): void;
}

export type { GridSearchProgressListener };
