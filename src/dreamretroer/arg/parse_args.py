import argparse
import os
import torch
import sys

parser = argparse.ArgumentParser()
#  gpu
parser.add_argument('--gpu', type=int, default=-1)


# random seed
parser.add_argument('--seed', type=int, default=1234)

# ================== one-step model ================ #
parser.add_argument('--mlp_model_dump',
                    default='dreamretroer/one_step_model/retro_star_value.ckpt')

parser.add_argument('--mlp_templates',
                    default='dreamretroer/one_step_model/template_rules.dat')
# ================================================== #
# search algs
parser.add_argument('--iterations', type=int, default=500)
parser.add_argument('--expansion_topk', type=int, default=50)
parser.add_argument('--viz', action='store_true')
parser.add_argument('--viz_dir', default='viz')

# collect or not
parser.add_argument('--collect_expe', default=True)
parser.add_argument('--experience_root', default='experience_dataset')
parser.add_argument('--experience_data', default='train_experience')


# eg_fn model
parser.add_argument('--mol_fp_dim', type=int, default=2048)
parser.add_argument('--template_fp_dim', type=int, default=2048)
parser.add_argument('--value_fp_dim', type=int, default=4096)


# one-step model
parser.add_argument('--fp_dim', type=int, default=2048)
parser.add_argument('--n_layers', type=int, default=1)
parser.add_argument('--latent_dim', type=int, default=256)

# train one round
parser.add_argument('--train_root', default='experience_dataset')
parser.add_argument('--train_data', default='train_experience')
parser.add_argument('--n_epochs', type=int, default=20)
parser.add_argument('--batch_size', type=int, default=512)
parser.add_argument('--lr', type=float, default=1e-3)
parser.add_argument('--save_epoch_int', type=int, default=1)
parser.add_argument('--save_folder', default='dreamretroer/one_step_model')

# evaluation
parser.add_argument('--use_value_fn', default=True)
parser.add_argument('--value_model', default='best_egn_for_emol_2.pt') 
parser.add_argument('--result_folder', default='dreamretroer/results')

# dataset
parser.add_argument('--test_mols', default='dreamretroer/dataset/retro190.txt') # 'retro190.txt'


parser.add_argument('--starting_molecules', default='dreamretroer/dataset/origin_dict.csv') # 'dataset/chembl.csv'
parser.add_argument('--reaction_data', default='dreamretroer/dataset/uspto.reactions.json')

# NOC
parser.add_argument('--noc_data', default='../dreamretroer/dataset/reaction_graph.gml')


args = parser.parse_args()

# setup device
os.environ["CUDA_VISIBLE_DEVICES"] = str(args.gpu)
