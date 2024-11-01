import json

import peptonizer

# Provided by Pyodide, required to send status updates from this thread to the main thread in JavaScript.
import js

class JSZeroLookaheadProgressListener(peptonizer.ZeroLookaheadProgressListener):
    def __init__(self, execution_id: int):
        self.execution_id = execution_id

    def graphs_updated(
        self,
        current_graph: int,
        total_graphs: int
    ):
        js.postMessage(json.dumps({
            "id": self.execution_id,
            "type": "progress",
            "progress_type": "graph",
            "current_graph": current_graph,
            "total_graphs": total_graphs
        }))

    def max_residual_updated(
        self,
        max_residual: float,
        tolerance: float,
    ):
        js.postMessage(json.dumps({
            "id": self.execution_id,
            "type": "progress",
            "progress_type": "max_residual",
            "max_residual": max_residual,
            "tolerance": tolerance
        }))

    def iterations_updated(
        self,
        current_iterations: int,
        total_iterations: int
    ):
        js.postMessage(json.dumps({
            "id": self.execution_id,
            "type": "progress",
            "progress_type": "iterations",
            "current_iterations": current_iterations,
            "total_iterations": total_iterations
        }))

print("Started running PepGM...")
graph = globals().get('graph')
alpha = globals().get('alpha')
beta = globals().get('beta')
prior = globals().get('prior')
execution_id = globals().get('execution_id')

pepgm_results = peptonizer.run_belief_propagation(
    graph,
    alpha,
    beta,
    True,
    prior,
    progress_listener=JSZeroLookaheadProgressListener(execution_id)
)

print("Successfully executed PepGM...")

# Now convert the results from PepGM into a list of taxon IDs and the corresponding score values.
final_scores = peptonizer.extract_taxon_scores(pepgm_results)

json.dumps(final_scores)