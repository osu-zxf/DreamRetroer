import pickle

routes = pickle.load(open('../dataset/routes_train.pkl', 'rb'))
# with open("../dataset/train_routes.txt", 'w') as f:
#     cnt = 0
#     for route in routes:
#         route_info = '**************** id=%d, len=%d ****************' % (cnt, len(route))
#         f.write(route_info)
#         f.write('\n')
#         for reaction in route:
#             f.write('%s\n' % reaction)
#         cnt += 1
    


length_counter = []
for i in range(100):
    length_counter.append(0)
for route in routes:
    if len(route) < 100:
        length_counter[len(route)]+=1
with open("../dataset/train_route_distribution.csv", 'w') as f:
    f.write("length,count\n")
    for i in range(len(length_counter)):
        length_info = '%d,%d' % (i,length_counter[i])
        f.write(length_info)
        f.write('\n')