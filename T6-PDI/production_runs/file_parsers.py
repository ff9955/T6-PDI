import numpy as np

def sh_log_parser(file_path):
    '''
    Parser made for reading the run-sh-1.log file, and filtering out the lines which show the index of the active state. It does this by only taking the lines whose final elements are integers, and passing them into an array of strings. We then get a 1D array of the active states at every timestep at which they were printed by X-SH, by just slicing the final column of the array of strings.
    '''

    with open(file_path) as full_file:

        filtered_list = []

        for single_line in full_file:

            split_line = single_line.split()

            try:
                active_state = int(split_line[-1])
                if len(split_line) == 4:
                    filtered_list.append(active_state)
                elif len(split_line) == 5:
                    filtered_list.append(active_state)
            except:
                continue

        active_states = np.array(filtered_list)
        return active_states

def xyz_parser(file_path, number_elements):
    '''
    reads .xyz or .txt file into an array of strings to be processed by another function, the line_elements argument refers to the number of elements your desired lines have, in case different lines in the file have different numebrs of elements
    '''

    with open(file_path) as full_file:

        filtered_list = []

        for single_line in full_file:

            split_line = single_line.split()

            if len(split_line) == number_elements:
                filtered_list.append(split_line)

    filtered_array = np.array(filtered_list)

    return filtered_array

def xyz_parse_section(file_path, number_elements, first_line, last_line):
    '''
    basically does what the above function does, but only between two line numbers that you can specify
    '''

    line_counter = 1
    filtered_list = []

    with open(file_path) as full_file:

        while line_counter <= last_line:

            single_line = full_file.readline()
            if line_counter >= first_line:

                split_line = single_line.split()
                if len(split_line) == number_elements:
                    filtered_list.append(split_line)

            line_counter += 1

    filtered_array = np.array(filtered_list)
    return filtered_array

def xyz_parse_first_section(file_path, number_elements):
    '''
    This parses the first section of a file that follows a certain pattern. For instance, if you tell the parser to only read in lines with 4 columns each, it will read in all lines with this pattern, UNTIL it reaches a line where it doesn't hold, then it will abort. However, it onyl does this after you've reached a line with the right pattern, meaning junk lines at the beginning of the file will not stop the parser.
    '''

    filtered_list = []

    with open(file_path) as full_file:

        single_line = True
        while single_line:

            single_line = full_file.readline()
            split_line = single_line.split()

            if len(split_line) == number_elements:
                filtered_list.append(split_line)
            elif filtered_list:
                break

    filtered_array = np.array(filtered_list)
    return filtered_array


