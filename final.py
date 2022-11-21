import re


def open_fasta_file(file_address):
    file = open(file_address, 'r')
    text = file.read()
    file.close()
    return text


def record_counter(file_address):
    txt = open_fasta_file(file_address)
    counter = txt.count('>')
    return counter


def dna_dict_creator(file_address):
    txt = open_fasta_file(file_address)
    dna_list = re.split('>', txt)
    dna_list = list(filter(None, dna_list))
    list_of_keys, new_dna_list = [], []
    totall_dna_in_fasta = {}

    for item in dna_list:
        dna = (re.split('\n', (item)))
        dna = list(filter(None, dna))
        list_of_keys.append(dna[0])
        dna = dna[1:]
        dnastr = ''
        dnastr = ''.join(dna)
        new_dna_list.append(dnastr)

    for i in range(0, len(list_of_keys)):
        totall_dna_in_fasta[list_of_keys[i][:30]] = new_dna_list[i]

    return totall_dna_in_fasta


def length_calculater(file_address):
    dna_dict = dna_dict_creator(file_address)
    dict_of_length = {}
    for item in dna_dict.keys():
        dict_of_length[item] = len(dna_dict[item])
    return dict_of_length


def longest_shortest(file_address):
    dna_length_dict = length_calculater(file_address)
    sorted_dna_length = sorted(dna_length_dict.items(), key=lambda x: x[1])
    longest = (sorted_dna_length[-1][0][:30], sorted_dna_length[-1][1])
    shortest = (sorted_dna_length[0][0][:30], sorted_dna_length[0][1])
    identical_shortest_len, identical_longest_len = [], []

    for i in range(0, len(sorted_dna_length)):
        longest_len = sorted_dna_length[-1][1]
        if sorted_dna_length[-i][1] == longest_len:
            identical_longest_len.append(sorted_dna_length[-i][0][:25])
        else:
            break

    for i in range(0, len(sorted_dna_length)):
        shortest_len = sorted_dna_length[0][1]
        if sorted_dna_length[i][1] == shortest_len:
            identical_shortest_len.append(sorted_dna_length[i][0][:25])
        else:
            break

    return shortest, longest, identical_shortest_len, identical_longest_len


def orf_finder(file_address, reading_frames):
    fasta_file = dna_dict_creator(file_address)
    dnas_keys = fasta_file.keys()

    dna_orf = {}
    for item in dnas_keys:
        dna = fasta_file[item]
        start_position = reading_frames - 1
        mark = 0
        start_index, stop_index, orf = [], [], []
        for i in range(start_position, len(dna), 3):
            if dna[i:i + 3] == 'ATG':
                start_index.append(i)
        for i in range(start_position, len(dna), 3):
            if dna[i:i + 3] in ["TAA", "TGA", "TAG"]:
                stop_index.append(i)
        for i in range(0, len(start_index)):
            for j in range(0, len(stop_index)):
                if start_index[i] < stop_index[j] and start_index[i] > mark:
                    orf.append(dna[start_index[i]:stop_index[j] + 3])
                    mark = stop_index[j] + 3
                    break
        dna_orf[item] = orf
    return dna_orf


def longest_orf_length(file_address, reading_frames):
    all_orfs = orf_finder(file_address, reading_frames)
    orf_keys = all_orfs.keys()
    orf_lengths = {}
    for identifier in orf_keys:
        orf = all_orfs[identifier]
        length = []
        for item in orf:
            length.append(len(item))
        length.sort(reverse=True)
        if len(length) > 0:
            orf_lengths[identifier] = length[0]
        else:
            orf_lengths[identifier] = 0
    return orf_lengths, max(orf_lengths.values())


def longest_orf_position(file_address, reading_frame):
    longest_orf_len = longest_orf_length(file_address, reading_frame)
    long_orfs = longest_orf_len[0]
    longest_orf_len = longest_orf_len[1]
    longest_orf_identifier = ''

    for item in long_orfs.keys():
        if long_orfs[item] == longest_orf_len:
            longest_orf_identifier = item

    all_dna = dna_dict_creator(file_address)
    all_orfs = orf_finder(file_address, reading_frame)

    longest_orf_in_fasta = all_orfs[longest_orf_identifier][0]
    for i in range(0, len(all_orfs[longest_orf_identifier])):
        if (len(all_orfs[longest_orf_identifier][i]) > len(longest_orf_in_fasta)):
            longest_orf_in_fasta = all_orfs[longest_orf_identifier][i]

    dna = all_dna[longest_orf_identifier]
    start_pos = dna.rfind(longest_orf_in_fasta)

    #
    # start_positions = {}
    #
    # for item in all_dna.keys():
    #     dna = all_dna[item]
    #     orfs = all_orfs[item]
    #     longest = ''
    #     for orf in orfs:
    #         if len(orf) > len(longest):
    #             longest = orf
    #     start_positions[item] = dna.rfind(longest)
    # return start_positions
    return start_pos + 1


def all_repeats(file_address, length):
    fasta_file = dna_dict_creator(file_address)
    dnas_keys = fasta_file.keys()
    repeats_dict = {}

    for item in dnas_keys:
        repeats_list, fragments = [], []
        dna = fasta_file[item]
        for i in range(0, len(dna)):
            fragments.append(dna[i:i + length])

        for piece in fragments:

            if len(piece) == length and dna.count(piece) > 1:
                repeats_list.append(dna.count(piece))
                repeats_dict[piece] = (dna.count(piece))

    return repeats_dict


def most_frequent_repeat(file_address, length):
    repeat_dict = all_repeats(file_address, length)
    keys = repeat_dict.keys()
    max_repeat = {}

    for key in keys:
        temp_list = []
        max_rep = max(repeat_dict[key])
        count = repeat_dict[key].count(max_rep)
        temp_list.append(max_rep)
        temp_list.append(count)
        max_repeat[key] = temp_list

    return max_repeat


# print(longest_orf_length('dna2.fasta', 3)[1])
# print(longest_orf_position('dna2.fasta',3))
# print(len(dna_dict_creator('dna.example.fasta')))
# print(longest_orf_position('dna.example.fasta',0))
# print(longest_orf_position('dna.example.fasta',0))
# print(most_frequent_repeat('dna.example.fasta',6))
# print(all_repeats('dna.example.fasta',6))

# orf = 'gi|142022655|gb|EQ086233.1|16'
# long_orfs = longest_orf_length('dna2.fasta',3)[0]
# for item in long_orfs.keys():
#     if orf in item:
#         print(long_orfs[item])


