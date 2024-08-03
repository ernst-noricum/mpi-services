import json
import uvicorn

from fastapi import FastAPI
from fastapi_healthcheck import HealthCheckFactory, healthCheckRoute

from dnachisel import *
from dnachisel.biotools import reverse_translate

from model import CodonOptimizeRequest, CodonOptimizeResponse


app = FastAPI()

# Add Health Checks
_healthChecks = HealthCheckFactory()
app.add_api_route('/health', endpoint=healthCheckRoute(factory=_healthChecks))


@app.get("/")
def get_root():
    return "Hello World"


@app.post("/codon_optimize", response_model=CodonOptimizeResponse)
async def codon_optimize(request: CodonOptimizeRequest):

    codon_optimized_seqs = list()

    for seq in request.sequences:

        dna_seq = seq.seq
        if seq.type == "Protein":
            # reverse translate to DNA first
            dna_seq = reverse_translate(seq.seq, randomize_codons=True, table='Standard')

        problem = DnaOptimizationProblem(
            sequence=dna_seq,
            constraints=[
                AvoidPattern("BsaI_site"),
                EnforceGCContent(mini=0.3, maxi=0.7, window=50),
                EnforceTranslation()
            ],
            objectives=[CodonOptimize(species='e_coli')]
        )

        # SOLVE THE CONSTRAINTS, OPTIMIZE WITH RESPECT TO THE OBJECTIVE
        problem.resolve_constraints()
        problem.optimize()

        # # PRINT SUMMARIES TO CHECK THAT CONSTRAINTS PASS
        # print(problem.constraints_text_summary())
        # print(problem.objectives_text_summary())

        # GET THE FINAL SEQUENCE (AS STRING OR ANNOTATED BIOPYTHON RECORDS)
        final_sequence = problem.sequence  # string
        # final_record = problem.to_record(with_sequence_edits=True)
        # print(final_record)

        codon_optimized_seqs.append({
            "original": {
                "name": seq.name,
                "type": seq.type,
                "seq": seq.seq
            },
            "codon_optimized": {
                "name": seq.name,
                "type": "DNA",
                "seq": final_sequence
            }
        })

    print(f"Response: {json.dumps(codon_optimized_seqs, indent=3)}")

    return {"sequences": codon_optimized_seqs}


if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=8000)