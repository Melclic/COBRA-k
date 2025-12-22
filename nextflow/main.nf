#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Build a COBRAK model with kinetic data and save as annotated SBML.
 * Container: cobrak:latest
 *
 * Example:
 *   nextflow run cobrak_enrich.nf \
 *     --input_model iML1515.xml \
 *     --database_data_folder database_data \
 *     --brenda_version 2025_1 \
 *     --base_species "Escherichia coli" \
 *     --prefer_brenda true \
 *     --do_model_fullsplit false \
 */

/* ----------------------------
 * Help message
 * ---------------------------- */
def helpMessage() {
  log.info """
╭────────────────────────────────────────────────────────────────────────╮
│                         COBRAK Model Enrichment                        │
╰────────────────────────────────────────────────────────────────────────╯

Usage:
  nextflow run cobrak_enrich.nf --input_model <SBML_FILE> [options] -with-docker

Required:
  --input_model <path>      Path to input SBML (.xml)

Options:
  --database_data_folder    Path to COBRAK database bundle (default: "database_data")
  --brenda_version          BRENDA version tag (default: "2025_1")
  --base_species            Species name (default: "Escherichia coli")
  --prefer_brenda           Prefer BRENDA data [true|false] (default: true)
  --do_model_fullsplit      Split model reactions [true|false] (default: false)
  --help                    Show this help

Example:
  nextflow run cobrak_enrich.nf --input_model iML1515.xml -with-docker
""".stripIndent()
}

/* ----------------------------
 * Parameters (with defaults)
 * ---------------------------- */
params.input_model            = null
params.base_species           = null
params.prefer_brenda          = true
params.do_model_fullsplit     = false
params.output_folder          = 'cobrak'
params.add_enzyme_constraints = false
params.combine_base_reactions = false
params.help                   = false

/* ----------------------------
 * Process
 * ---------------------------- */
process build_cobrak_model {
  errorStrategy 'terminate'

  container "cobrak:latest"

  input:
    path model_file
    val  base_species
    val  prefer_brenda
    val  do_fullsplit

  output:
    path "cobrak_model.json", emit: cobrak_model

  script:
  """
#!/usr/bin/env python3
from cobrak.model_instantiation import get_cobrak_model_with_kinetic_data_from_sbml_model_alone
from cobrak.io import save_cobrak_model_as_annotated_sbml_model
from cobra.io import read_sbml_model

model_path = "${model_file}"
prefer_brenda = True if "${prefer_brenda}".lower() in ("true", "True", "TRUE") else False
do_model_fullsplit = True if "${do_fullsplit}".lower() in ("true", "True", "TRUE") else False
base_species = None if str("${base_species}") in ("None", "", "null", "NaN", "nan") else str("${base_species}")

add_enzyme_constraints = True if "${params.add_enzyme_constraints}".lower() in ("true", "True", "TRUE") else False
combine_base_reactions = True if "${params.combine_base_reactions}".lower() in ("true", "True", "TRUE") else False

#try to recover the base_species name or taxonomy id
if not base_species:
    tmp_model = read_sbml_model(model_path) 
    base_species = tmp_model.annotation.get('taxonomy')
    tmp_model = None

m = get_cobrak_model_with_kinetic_data_from_sbml_model_alone(
    sbml_path=model_path,
    database_data_folder="database_data",
    brenda_version="2025_1",
    base_species=base_species,
    prefer_brenda=prefer_brenda,
    do_model_fullsplit=do_model_fullsplit,
)

from cobrak.io import json_write
json_write(
    path="cobrak_enriched_model.json",
    json_data=m,
)
  """
}

process run_cobrak_lp {
  errorStrategy 'terminate'
  container "cobrak:latest"

  input:
    path model_file
    val  mode
    val  objective_str
    val  objective_json
    val  sense
    val  output_folder
    val  with_flux_sum_var
    val  with_enzyme_constraints

  output:
    path "result.json", optional: true, emit: result_json
    path "variability.json", optional: true, emit: variability_json

  script:
  """
#!/usr/bin/env python3
import json
from pathlib import Path
from typing import Any, Dict, Union, Optional

from cobrak.io import json_load, json_write
from cobrak.dataclasses import Model
from cobrak.lps import perform_lp_optimization, perform_lp_variability_analysis
from cobrak.constants import FLUX_SUM_VAR_ID

model_path = Path("${model_file}").resolve()

mode = "${mode}".strip().lower()
sense = int("${sense}")

with_flux_sum_var = str("${with_flux_sum_var}").lower() in {"1","true","yes","y"}
with_enzyme_constraints = str("${with_enzyme_constraints}").lower() in {"1","true","yes","y"}

# objective can come in as either:
# - objective_str: a single reaction id (string), OR empty
# - objective_json: a JSON dict string like {"Glycolysis": -1.5, "Overflow": 2.0}, OR empty
objective: Optional[Union[str, Dict[str, float]]] = None

objective_str = "${objective_str}".strip()
objective_json = "${objective_json}".strip()

if objective_json:
    objective = json.loads(objective_json)
elif objective_str:
    objective = objective_str

cobrak_model: Model = json_load(
    path=str(model_path),
    dataclass_type=Model,
)

result_path = "result.json"
variability_path = "variability.json"

if mode == "fba":
    if objective is None:
        raise ValueError("objective is required for mode=fba")
    res = perform_lp_optimization(
        cobrak_model,
        objective,
        sense,
    )
    json_write(str(result_path), res)

elif mode == "pfba":
    if objective is None:
        raise ValueError("objective is required for mode=pfba")

    # 1) optimize objective first
    fba_res = perform_lp_optimization(
        cobrak_model,
        objective,
        sense,
    )

    # choose a reaction id to pin as min flux
    if isinstance(objective, str):
        obj_rxn = objective
    else:
        obj_rxn = next(iter(objective.keys()))

    cobrak_model.reactions[obj_rxn].min_flux = float(fba_res[obj_rxn])

    # 2) minimize flux sum
    pfba_res = perform_lp_optimization(
        cobrak_model,
        FLUX_SUM_VAR_ID,
        -1,
        with_flux_sum_var=True,
    )

    # reset
    cobrak_model.reactions[obj_rxn].min_flux = 0.0

    json_write(str(result_path), {"fba": fba_res, "pfba": pfba_res})

elif mode == "fva":
    var_res = perform_lp_variability_analysis(
        cobrak_model,
    )
    json_write(str(variability_path), var_res)

elif mode == "ecfba":
    if objective is None:
        raise ValueError("objective is required for mode=ecfba")
    res = perform_lp_optimization(
        cobrak_model,
        objective,
        sense,
        with_enzyme_constraints=True,
    )
    json_write(str(result_path), res)

elif mode == "ecfva":
    var_res = perform_lp_variability_analysis(
        cobrak_model,
    )
    json_write(str(variability_path), var_res)

else:
    raise ValueError(f"Unknown mode: {mode}")
  """
}

/* ----------------------------
 * Workflow
 * ---------------------------- */
workflow {
    if (params.help || !params.input_model) {
        helpMessage()
        exit 0
    }

    ch_model = channel.fromPath(params.input_model)
    ch_species = channel.value(params.base_species)
    ch_pref = channel.value(params.prefer_brenda.toString())
    ch_split = channel.value(params.do_model_fullsplit.toString())
    ch_out = channel.value(params.output_folder)

    build_cobrak_model(
        ch_model,
        ch_species,
        ch_pref,
        ch_split,
        ch_out
    )
}

