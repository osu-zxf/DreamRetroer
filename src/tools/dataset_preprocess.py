import pandas as pd

all_data_file = "../dataset/uspto_reactions.csv"
df = pd.read_csv(all_data_file)

train_set = df.sample(frac=0.8, random_state=60)
train_set.to_csv("../dataset/uspto_reactions_train.csv", index=False)
eval_set = df.sample(frac=0.1, random_state=80)
eval_set.to_csv("../dataset/uspto_reactions_eval.csv", index=False)
test_set = df.sample(frac=0.1, random_state=100)
test_set.to_csv("../dataset/uspto_reactions_test.csv", index=False)