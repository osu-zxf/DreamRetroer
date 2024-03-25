common_path_list = []
with open('../dataset/inter_list.txt', 'r') as f:
    cnt = 0
    while True:
        cnt += 1
        # Get next line from file
        line = f.readline().strip()
        # if line is empty
        # or end of file is reached
        if not line:
            break
        if cnt % 2 == 0:
            path_str_list = (line[1:-1]).split(',') # [xxxxx]
            for path_str in path_str_list:
                path = (path_str.strip())[1:-1]  # 'xxxxx'
                if len(path.split('->')) == 2:
                    common_path_list.append(path)

common_path_list = list(set(common_path_list))

with open('../dataset/common_path_chain_2.txt', 'w') as f:
    for path in common_path_list:
        f.write('%s\n' % path)