import os
import tarfile
import tempfile
import logging
import sys


from cobrak.model_instantiation import get_cobrak_model_with_kinetic_data_from_sbml_model_alone
from cobrak.io import save_cobrak_model_as_annotated_sbml_model

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

bigg_model_species = {
    "iAF1260.xml": "Escherichia coli",
    "iAF987.xml": "Escherichia coli",
    "iJN1463.xml": "Escherichia coli",
    "iJN746.xml": "Synechocystis sp",
    "iJO1366.xml": "Escherichia coli",
    "iJR904.xml": "Escherichia coli",
    "iML1515.xml": "Escherichia coli",
    "iMM904.xml": "Saccharomyces cerevisiae",
    "iND750.xml": "Saccharomyces cerevisiae",
    "iYO844.xml": "Yersinia pestis",
}


def process_models_from_tar(
    tar_path: str,
    database_data_folder: str = "/COBRA-k/database_data",
    brenda_version: str = "2025_1",
    prefer_brenda: bool = True,
    do_model_fullsplit: bool = False,
    output_dir: str = "processed_models"
):
    """Extract all SBML models from a .tar.gz file and process them with COBRA-k.

    Args:
        tar_path (str): Path to models.tar.gz file.
        database_data_folder (str): Folder containing COBRA-k database data.
        brenda_version (str): Version of the BRENDA database (e.g., '2025_1').
        base_species (str): Base species name for kinetic data.
        prefer_brenda (bool): Prefer BRENDA database entries if True.
        do_model_fullsplit (bool): Whether to fully split the model reactions.
        output_dir (str): Directory to save processed COBRA-k annotated models.
    """
    os.makedirs(output_dir, exist_ok=True)

    with tempfile.TemporaryDirectory() as tmpdir:
        logging.info(f"Extracting {tar_path} → {tmpdir}")
        with tarfile.open(tar_path, "r:gz") as tar:
            tar.extractall(tmpdir)

        for root, _, files in os.walk(tmpdir):
            for filename in files:
                if filename.endswith(".xml"):
                    sbml_path = os.path.join(root, filename)
                    logging.info(f"Processing model: {filename}")
                    print(f"Processing model: {filename}", flush=True)

                    try:
                        cobrak_model = get_cobrak_model_with_kinetic_data_from_sbml_model_alone(
                            sbml_path=sbml_path,
                            database_data_folder=database_data_folder,
                            brenda_version=brenda_version,
                            base_species=bigg_model_species.get(os.path.basename(sbml_path)),
                            prefer_brenda=prefer_brenda,
                            do_model_fullsplit=do_model_fullsplit,
                        )

                        output_path = os.path.join(output_dir, f"cobrak_{filename}")
                        save_cobrak_model_as_annotated_sbml_model(
                            cobrak_model,
                            filepath=output_path,
                            combine_base_reactions=False,
                            add_enzyme_constraints=False,
                        )
                        logging.info(f"Saved processed model to: {output_path}")
                        print(f"Saved processed model to: {output_path}", flush=True)

                    except Exception as e:
                        logging.error(f"Failed to process {filename}: {e}")
                        print(f"Failed to process {filename}: {e}", flush=True)

if __name__ == "__main__":
    process_models_from_tar("models.tar.gz", output_dir='/COBRA-k/enriched_models/')
