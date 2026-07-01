
# %pip install pykeen
# %pip install torch torchvision

import random
import numpy as np
import torch
import pandas as pd
import numpy as np
from pykeen.triples import TriplesFactory
from pykeen.hpo import hpo_pipeline
from pykeen import predict
from pykeen.pipeline import pipeline_from_config
from pykeen.utils import set_random_seed
import logging 
import os

randomseed = 48 # Change this for each run to get different results
torch.manual_seed(randomseed)
random.seed(randomseed)
np.random.seed(randomseed)
set_random_seed(randomseed)

folder_path = f'seed{randomseed}/alzkb_rotate_with_br/'
os.makedirs(folder_path, exist_ok=True)

data = pd.read_pickle('/home/khanjanin/miniconda3/envs/neda_myenv_PyKEEN_AlzKB_ABA/alzkb_br/alzkb_br_triples_memgraph_v2.pkl')

# Configure logging to save to a file
log_file = folder_path+'hpo_pipeline.log'

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s:%(message)s', handlers=[
    logging.FileHandler(log_file),
    logging.StreamHandler()
])
logger = logging.getLogger(__name__)

try:
    logger.info("Starting the process")
    logger.info(data.shape)
    logger.info(set(data['type']))
    
    # Generating triples from graph
    if 'source' not in data.columns or 'type' not in data.columns or 'target' not in data.columns:
        logger.error("Data does not have required columns: 'source', 'type', 'target'")
        raise ValueError("Data must have 'source', 'type', and 'target' columns")

    if data.isnull().any().any():
        logger.warning("Data contains null values. These will be dropped.")
        data = data.dropna()

    logger.info("Creating TriplesFactory from labeled triples")
    tf = TriplesFactory.from_labeled_triples(triples = data[["source", "type", "target"]].values)
    # tf.to_path_binary(folder_path+'triples')
    logger.info("Created TriplesFactory from labeled triples with %d entities and %d relations", tf.num_entities, tf.num_relations)

    # Random split
    training, testing, validation = tf.split([.8, .1, .1], random_state=0)
    logger.info(f"Split data into training, testing, and validation sets: "
                f"Training: {training.num_triples}, Testing: {testing.num_triples}, Validation: {validation.num_triples}")

    # Training RotatE model
    logger.info("Starting training of RotatE model")
    result_rotate = hpo_pipeline(
        training=training,
        testing=testing,
        validation=validation,
        model='RotatE',
        model_kwargs=dict(
            random_seed=42
        ),
        stopper='early',
        epochs=100,
        device='gpu',
        evaluator='RankBasedEvaluator',
        n_trials=100,
        save_model_directory=folder_path+'rotate_models/',
    ) #device can be cpu or gpu
    logger.info("Finished RotatE pipeline")

    # save the model
    result_rotate.save_to_directory(folder_path+'alzkb_result_rotate_hpo')
    logger.info("Saved the model to directory")

    # get best config from hyperparameter optimization
    best_cfg = result_rotate._get_best_study_config()

    best_cfg['pipeline']['training'] = training
    best_cfg['pipeline']['testing'] = testing
    best_cfg['pipeline']['validation'] = validation
    
    rotate_pipeline = pipeline_from_config(best_cfg)
    logger.info("Finished training RotatE model")

    # save the model
    rotate_pipeline.save_to_directory(folder_path+'alzkb_result_rotate')
    logger.info("Saved the model to directory")

    rank_metrics = rotate_pipeline.metric_results.to_df()
    
    # Extracting the values
    hits_at_1 = float(rank_metrics[rank_metrics['Metric'] == 'hits_at_1'].head(1)['Value'])
    hits_at_3 = float(rank_metrics[rank_metrics['Metric'] == 'hits_at_3'].head(1)['Value'])
    hits_at_10 = float(rank_metrics[rank_metrics['Metric'] == 'hits_at_10'].head(1)['Value'])
    adjusted_inverse_harmonic_mean_rank = float(rank_metrics[rank_metrics['Metric'] == 'adjusted_inverse_harmonic_mean_rank'].head(1)['Value'])

    # Formatting the values into a string
    results = f"{hits_at_1} {hits_at_3} {hits_at_10} {adjusted_inverse_harmonic_mean_rank}"
    logger.info(f"Extracted rank metrics: {results}")

    # Saving the results to a text file
    with open(folder_path+'rotate_results.txt', "w") as file:
        file.write(results)
    logger.info("Saved results to text file")

    rank_metrics.to_csv(folder_path+'rank_metric_rotate.csv')
    logger.info("Saved rank metrics to CSV file")

    node_id = '233508'  #node id for Alzheimer's Disease node

    drugs_df = predict.predict_target(model = rotate_pipeline.model, relation = 'DRUGTREATSDISEASE', tail = node_id, triples_factory=training).df
    drugs_df.to_csv(folder_path+'rotate_head_drug_prediction.csv', index=False)
    logger.info("Saved drug predictions to CSV file")

    # Extract entity and relation embeddings
    entity_embeddings = rotate_pipeline.model.entity_representations[0](indices=None).detach().cpu().numpy()
    relation_embeddings = rotate_pipeline.model.relation_representations[0](indices=None).detach().cpu().numpy()
    logger.info("Extracted entity and relation embeddings")

    # Save embeddings to CSV
    entity_df = pd.DataFrame(entity_embeddings)
    entity_df.to_csv(folder_path+'rotate_entity_embeddings.csv', index=False)
    relation_df = pd.DataFrame(relation_embeddings)
    relation_df.to_csv(folder_path+'rotate_relation_embeddings.csv', index=False)
    logger.info("Saved embeddings to CSV files")

    # Optional: Save embeddings to NumPy files
    np.save(folder_path+'rotate_entity_embeddings.npy', entity_embeddings)
    np.save(folder_path+'rotate_relation_embeddings.npy', relation_embeddings)
    logger.info("Saved embeddings to NumPy files")

    logger.info("Process completed successfully")


except Exception as e:
    logger.error(f"An error occurred: {e}", exc_info=True)
