import numpy as np

HOMO_wf = []

new_FMO_file = open('hole-density.cube', 'a')

def round_4dp(number):

    rounded_number = f'{number:.4e}'
    return float(rounded_number)

vector_round = np.vectorize(round_4dp)

with open('hole_crystal_286.cube') as LUMO_cube:
    with open('hole_crystal_286.cube') as HOMO_cube:
       
        for j in range(50):
            HOMO_cube.readline()
            LUMO_cube.readline() 

        for line in HOMO_cube:
            split_line = line.split()
            split_line = list(map(lambda x: float(x), split_line))
            
            HOMO_wf.append(split_line)

        LUMO_lines = LUMO_cube.readlines()

        for j in range(len(LUMO_lines)):
            LUMO_line = LUMO_lines[j]
            split_line = LUMO_line.split()
            split_line = list(map(lambda x: float(x), split_line))

            if not split_line:
                new_FMO_file.write('\n')
                continue

            LUMO_array = np.array(split_line)
            HOMO_array = np.array(HOMO_wf[j])

            conjugated_array = LUMO_array*HOMO_array
            conjugated_array = vector_round(conjugated_array)

            conjugated_string = str(conjugated_array)[1:-1]

            new_FMO_file.write(conjugated_string + '\n')

new_FMO_file.close()
