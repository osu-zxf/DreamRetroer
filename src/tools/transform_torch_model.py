import torch
print(torch.__version__)
model = torch.load('./saved_EG_fn/best_egn_for_emol.pt')
# print(model)
torch.save(model, './saved_EG_fn/best_egn_for_emol_2.pt', _use_new_zipfile_serialization=False)