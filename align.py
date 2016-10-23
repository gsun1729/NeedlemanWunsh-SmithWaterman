__author__ = 'Gordon Sun'
# for BIOE 214 Fall 2015
import math
import sys
# global variables
global all_aligned_seq, gap_char, seq_1, seq_2, align_mode, score, M, Ix, Iy, M_dependency, Ix_dependency, Iy_dependency
all_aligned_seq = []  # all of the aligned sequences will be printed in here
gap_char = "_"

# input_filename = "alignment3.input.txt"  # for testing purposes only
# output_filename = "outout.txt"
# reads in input file name and output file name

input_filename = sys.argv[-2]
output_filename = sys.argv[-1]


# function generates a 2D array full of zeros
# I: (y,x) dimensions in a tuple
# O: 2D array of zero
def zeros(dimensions):
    zeros_array = []
    for x in range(dimensions[0]):
        zeros_array.append([])
        for y in range(dimensions[1]):
            zeros_array[-1].append(0.0)
    return zeros_array


# function generates a 2D array full of "000", for dependency tracing
# I: (y,x) dimensions in a tuple
# O: 2D array of "000"
def zeros_dependency(dimensions):
    zeros_array = []
    for x in range(dimensions[0]):
        zeros_array.append([])
        for y in range(dimensions[1]):
            zeros_array[-1].append("000")
    return zeros_array


# function to print a list in 2D format
# I: a, a is a 2D list
# O: printed list in 2D format.
def print_2d_list(a):
    if a == []:
        # So we don't crash accessing a[0]
        print []
        return
    rows, cols = len(a), len(a[0])
    field_width = max_item_length(a)
    print "[ ",
    for row in xrange(rows):
        if row > 0: print "\n  ",
        print "[ ",
        for col in xrange(cols):
            if col > 0: print ",",
            # The next 2 lines print a[row][col] with the given field_width
            formatted_array = "%" + str(field_width) + "s"
            print formatted_array % str(a[row][col]),
        print "]",
    print "]"


# function to help 2d list print; determines field width
# I:  array to be printed (a)
# O:  field width or longest item length in the array
def max_item_length(a):
    max_len, rows, cols = 0, len(a), len(a[0])
    for row in xrange(rows):
        for col in xrange(cols):
            max_len = max(max_len, len(str(a[row][col])))
    return max_len


# Defining a function that will return the score of two matching letters ONLY,
# this is used to determine the score only if matched (gaps not accounted for in this)
# Accepts two letters and the relevant score matrix as input and outputs score
# I: array score; char letter1; char letter2
# O: float score
def return_match_score(letter1, letter2):
    for j in xrange(len(score)):
        if letter1 == score[j][0] and letter2 == score[j][1]:
            return float(score[j][len(score[0]) - 1])


# Given an index, returns the relevant letter with the associated index in the input sequence
# Starts at 0; i.e. ABCD if letter index is 0 then returns A
# I: k is index (starting from 0) and seq_num is the sequence identifier
# O: single letter (char) of that sequence
def extract_seq_letter(k, seq_num):
    # Sequence 1/A is defined as seq_num=0 and
    # Sequence 2/B is defined as seq_num=1
    # k is defined as the index of the letter in the sequence
    k -= 1
    if k >= 0:
        if seq_num == 0:
            return seq_1[k]
        else:
            return seq_2[k]
    else:
        return " "


# Defining M_i_j function
# i-->Seq_1/A
# j-->Seq_2/B
# reads from M, Ix, Iy arrays to fill next cells
# Returns the max score as well as the array it was derived from
# I: i is index of sequence1, j is index of sequence2
# O: populates F score dependency matrix with dependencies, as demonstrated below:
# M=index 1
# Ix=index2
# Iy=index3
# product of combinations of 3 values indicates status of dependency, so if it returns 111, M, Iy,Ix are all dependencies,
# while if it returns 011, Ix, Iy are dependencies and so forth
def M_i_j(i, j):
    S = return_match_score(seq_1[i - 1], seq_2[j - 1])  # -1 because of M offset from corner by one row and column
    A = float(M[i - 1][j - 1] + S)
    B = float(Ix[i - 1][j - 1] + S)
    C = float(Iy[i - 1][j - 1] + S)
    if align_mode == 1:
        if A < 0:
            A = 0.0
        if B < 0:
            B = 0.0
        if C < 0:
            C = 0.0
    max_val = round(max(A, B, C), 1)

    if round(max(A, B, C), 1) == round(A, 1):
        if round(A, 1) == round(B, 1) and round(A, 1) == round(C, 1):
            return max_val, "111"
        elif round(A, 1) == round(C, 1):
            return max_val, "101"
        elif round(A, 1) == round(B, 1):
            return max_val, "110"
        else:
            return max_val, "100"
    elif round(max(A, B, C), 1) == round(B, 1):
        if round(B, 1) == round(A, 1) and round(B, 1) == round(C, 1):
            return max_val, "111"
        elif round(B, 1) == round(A, 1):
            return max_val, "110"  # "M", "Ix"
        elif round(B, 1) == round(C, 1):
            return max_val, "011"  # "Ix", "Iy"
        else:
            return max_val, "010"  # "Ix"
    else:
        if round(C, 5) == round(A, 5) and round(C, 5) == round(B, 5):
            return max_val, "111"  # M,Ix,Iy
        elif round(C, 5) == round(A, 5):
            return max_val, "101"  # "M", "Iy")
        elif round(C, 5) == round(B, 5):
            return max_val, "011"  # "Ix", "Iy")
        else:
            return max_val, "001"  # "Iy")


# Defining I_x function, reads from M and Ix arrays to populate the next cells
# Returns the max score as well as the array it was derived from
# I: i is index of sequence1, j is index of sequence2
# O: populates F score dependency matrix with dependencies, as demonstrated below:
# M=index 1
# Ix=index2
# Iy=index3
# product of combinations of 3 values indicates status of dependency, so if it returns 111, M, Iy,Ix are all dependencies,
# while if it returns 011, Ix, Iy are dependencies and so forth
def I_x(i, j):
    A = float(M[i - 1][j] - seq2_dx)
    B = float(Ix[i - 1][j] - seq2_ex)
    if align_mode == 1:
        if A < 0:
            A = 0.0
        if B < 0:
            B = 0.0
    max_val = round(max(A, B), 1)
    if round(max(A, B), 1) == round(A, 1):
        if round(A, 1) == round(B, 1):
            return max_val, "110"  # "M", "Ix")
        else:
            return max_val, "100"  # "M")
    else:
        if round(A, 1) == round(B, 1):
            return max_val, "110"  # "M", "Ix")
        else:
            return max_val, "010"  # "Ix")


# Defining I_y function, reads from M and Iy arrays to populate the next cells
# Returns the max score as well as the array it was derived from
# I: i is index of sequence1, j is index of sequence2
# O: populates F score dependency matrix with dependencies, as demonstrated below:
# M=index 1
# Ix=index2
# Iy=index3
# product of combinations of 3 values indicates status of dependency, so if it returns 111, M, Iy,Ix are all dependencies,
# while if it returns 011, Ix, Iy are dependencies and so forth
def I_y(i, j):
    A = float(M[i][j - 1] - seq1_dx)
    B = float(Iy[i][j - 1] - seq1_ex)
    if align_mode == 1:
        if A < 0:
            A = 0.0
        if B < 0:
            B = 0.0
    max_val = round(max(A, B), 1)
    if round(max(A, B), 1) == round(A, 1):
        if round(A, 1) == round(B, 1):
            return max_val, "101"  # "M", "Iy")
        else:
            return max_val, "100"  # "M")
    else:
        if round(A, 1) == round(B, 1):
            return max_val, "101"  # "M", "Iy")
        else:
            return max_val, "001"  # "Iy")


# Returns the maximum score found in the array depending on the alignment mode
# searches through M, Ix, Iy
# I: none; reads the alignment mode based on global variable (1 for local, 0 for global
# O: returns the max score out of all three arrays (M, Ix, Iy)
def find_max_score():
    max_score = 0
    if align_mode == 1:
        for x3 in xrange(1, len(M)):
            for y3 in xrange(1, len(M[0])):
                if max(M[x3][y3], Ix[x3][y3], Iy[x3][y3]) > max_score:
                    max_score = max(M[x3][y3], Ix[x3][y3], Iy[x3][y3])
    else:
        for botrow in xrange(1, len(M[0])):
            if max(M[len(M) - 1][botrow], Iy[len(Iy) - 1][botrow], Ix[len(Ix) - 1][botrow]) > max_score:
                max_score = max(M[len(M) - 1][botrow], Iy[len(Iy) - 1][botrow], Ix[len(Ix) - 1][botrow])
        for lastcol in xrange(1, len(M)):
            if max(M[lastcol][len(M[0]) - 1], Iy[lastcol][len(Iy[0]) - 1], Ix[lastcol][len(Ix[0]) - 1]) > max_score:
                max_score = max(M[lastcol][len(M[0]) - 1], Iy[lastcol][len(Iy[0]) - 1], Ix[lastcol][len(Ix[0]) - 1])
    return max_score


# Takes a score and determines the indexes at which the score shows up
# outputs a tuple (array name, I, J) with the array name and the index at which the score was found
# if there are multiple scores are the max but have different indices, they are arranged in an output 2D list of the form
# (array,i1,j1),(array,i2,j2)...
# Returns the array in which the score is found as well
# I: float score value (score)
# O: array of matrix identities and indexes where the score has been found in the format of
# (array,i1,j1),(array,i2,j2)...
def find_max_score_index(score):
    indexes = []
    if align_mode == 1:
        for x6 in xrange(1, len(M)):
            for y6 in xrange(1, len(M[0])):
                if M[x6][y6] == score:
                    indexes.append(("M", x6, y6))
                if Ix[x6][y6] == score:
                    indexes.append(("Ix", x6, y6))
                if Iy[x6][y6] == score:
                    indexes.append(("Iy", x6, y6))
    else:
        for d1 in xrange(1, len(M[0])):
            if M[len(M) - 1][d1] == score:
                indexes.append(("M", len(M) - 1, d1))
            if Ix[len(Ix) - 1][d1] == score:
                indexes.append(("Ix", len(Ix) - 1, d1))
            if Iy[len(Iy) - 1][d1] == score:
                indexes.append(("Iy", len(Iy) - 1, d1))
        for e1 in xrange(1, len(M) - 1):  # this one has a -1 to prevent double counting the bottom right cell
            if M[e1][len(M[0]) - 1] == score:
                indexes.append(("M", e1, len(M[0]) - 1))
            if Ix[e1][len(Ix[0]) - 1] == score:
                indexes.append(("Ix", e1, len(Ix[0]) - 1))
            if Iy[e1][len(Iy[0]) - 1] == score:
                indexes.append(("Iy", e1, len(Iy[0]) - 1))
    return indexes


# Recursive traceback given a starting point determines the next cell and sequence alignment until it reaches (0,0)
# I: array to start in in string format ("M", "Ix", or "Iy"), index of where to start in the array (i and j),
# and the sequence as far as it has been determined in the form of a tuple (seq1, seq2)
# O: returns nothing, just determines the sequences based upon the populated score arrays M Ix and Iy
def traceback(start_array, i, j, seqs):
    s1, s2 = seqs[0], seqs[1]
    if align_mode == 0:  # GLOBAL
        if start_array == "M":
            temp = M_dependency[i][j]
            s1 = extract_seq_letter(i, 0) + s1
            s2 = extract_seq_letter(j, 1) + s2
            if i == 0 or j == 0:
                all_aligned_seq.append((s1, s2))
                return
            else:
                if temp[0] == "1":
                    traceback("M", i - 1, j - 1, (s1, s2))
                if temp[1] == "1":
                    traceback("Ix", i - 1, j - 1, (s1, s2))
                if temp[2] == "1":
                    traceback("Iy", i - 1, j - 1, (s1, s2))

        elif start_array == "Ix":
            temp = Ix_dependency[i][j]
            s1 = extract_seq_letter(i, 0) + s1
            s2 = gap_char + s2
            if temp[0] == "1":
                traceback("M", i - 1, j, (s1, s2))
            if temp[1] == "1":
                traceback("Ix", i - 1, j, (s1, s2))
        elif start_array == "Iy":
            temp = Iy_dependency[i][j]
            s1 = gap_char + s1
            s2 = extract_seq_letter(j, 1) + s2
            if temp[0] == "1":
                traceback("M", i, j - 1, (s1, s2))
            if temp[2] == "1":
                traceback("Iy", i, j - 1, (s1, s2))
    else:  # local
        if start_array == "M":
            temp = M_dependency[i][j]
            s1 = extract_seq_letter(i, 0) + s1
            s2 = extract_seq_letter(j, 1) + s2
            if (temp == "110" and M[i - 1][j - 1] == 0 and Ix[i - 1][j - 1] == 0) or (
                                temp == "101" and M[i - 1][j - 1] == 0 and Iy[i - 1][j - 1] == 0) or (
                                    temp == "111" and M[i - 1][j - 1] == 0 and Iy[i - 1][j - 1] == 0 and Ix[i - 1][
                            j - 1] == 0):
                all_aligned_seq.append((s1, s2))
            else:
                if temp[0] == "1":
                    traceback("M", i - 1, j - 1, (s1, s2))
                if temp[1] == "1":
                    traceback("Ix", i - 1, j - 1, (s1, s2))
                if temp[2] == "1":
                    traceback("Iy", i - 1, j - 1, (s1, s2))

        elif start_array == "Ix":
            temp = Ix_dependency[i][j]
            s1 = extract_seq_letter(i, 0) + s1
            s2 = gap_char + s2
            # print s1, s2
            if (temp == "110" and M[i - 1][j] == 0 and Ix[i - 1][j] == 0):
                traceback("M", i - 1, j, (s1, s2))
            else:
                if temp[0] == "1":
                    traceback("M", i - 1, j, (s1, s2))
                if temp[1] == "1":
                    traceback("Ix", i - 1, j, (s1, s2))
        elif start_array == "Iy":
            temp = Iy_dependency[i][j]
            s1 = gap_char + s1
            s2 = extract_seq_letter(j, 1) + s2
            # print s1, s2
            if (temp == "101" and M[i][j - 1] == 0 and Iy[i][j - 1] == 0):
                traceback("M", i, j - 1, (s1, s2))
            else:
                if temp[0] == "1":
                    traceback("M", i, j - 1, (s1, s2))
                if temp[2] == "1":
                    traceback("Iy", i, j - 1, (s1, s2))


# ===========================================================================================
# File import mechanism, accepts on the command line two arguments, one for the input file and the other for the
# output file. Syntax of command MUST be in the order of INPUT file then OUTPUT file name
# EX: align.py inputfile.txt outputfile.txt
# note files must be in txt format

f = open(input_filename)
# Data is an array containing every line in the input file.
data = f.readlines()
# Extract sequences from the data file and strip the newline character from the list.
# produces each sequence into a list format
seq_1 = list(data[0].rstrip('\n'))
seq_2 = list(data[1].rstrip('\n'))
align_mode = int(data[2].rstrip('\n'))
score = []
# 1=local,0=global
# Extract penalties and put into variables
penalties = data[3].rstrip('\n')
list_pen = penalties.split()
seq1_dx = float(list_pen[0])
seq1_ex = float(list_pen[1])
seq2_dx = float(list_pen[2])
seq2_ex = float(list_pen[3])
# Extract alphabet lengths , characters we don't need as they are provided in the score array
seq1_alpha_len = len(data[5].rstrip('\n'))
seq2_alpha_len = len(data[7].rstrip('\n'))
# Extract the score array into three column format [letter1, letter2, score]
for i in xrange(8, len(data)):
    string_n = data[i].rstrip('\n').split()  # c is a temporary variable
    score.append(string_n[2:len(string_n)])

# Make 3 blank score arrays
# Sequence A axis aligned along the LEFT
# Sequence B axis aligned along TOP
# Hence Iy will be horizontal arrow
# Ix will be vertical arrow
# M will be diagonal arrow

# Note the top row and first column of M, Ix, and Iy will all initialize and remain 0
# create blank arrays for populating
M = zeros((len(seq_1) + 1, len(seq_2) + 1))
Ix = zeros((len(seq_1) + 1, len(seq_2) + 1))
Iy = zeros((len(seq_1) + 1, len(seq_2) + 1))
M_dependency = zeros_dependency((len(seq_1) + 1, len(seq_2) + 1))
Ix_dependency = zeros_dependency((len(seq_1) + 1, len(seq_2) + 1))
Iy_dependency = zeros_dependency((len(seq_1) + 1, len(seq_2) + 1))

for x1 in xrange(1, len(M)):  # col travel
    for y1 in xrange(1, len(M[0])):  # row travel
        # F G and H consist of tuple (score, 000) <dependencies
        F = M_i_j(x1, y1)
        G = I_x(x1, y1)
        H = I_y(x1, y1)

        M[x1][y1] = F[0]
        Ix[x1][y1] = G[0]
        Iy[x1][y1] = H[0]

        M_dependency[x1][y1] = str(F[1])
        Ix_dependency[x1][y1] = str(G[1])
        Iy_dependency[x1][y1] = str(H[1])

scores_found = find_max_score()
array_of_max_scores = find_max_score_index(scores_found)
traceback(array_of_max_scores[0][0], array_of_max_scores[0][1], array_of_max_scores[0][2], ("", ""))

numb_align, pair = len(all_aligned_seq), len(all_aligned_seq[0])
# file output, data is written to a txt file with a user determined name
# t scores_found
# print_2d_list(all_aligned_seq)
# print "========="
# print_2d_list(M)
# print "========="
# print_2d_list(Ix)
outputfile = open(output_filename, "w")
outputfile.write(str(scores_found) + "\n\n")
for printx in xrange(0, numb_align):
    for printy in xrange(0, pair):
        outputfile.write(str(all_aligned_seq[printx][printy]) + "\n")
    outputfile.write("\n")
outputfile.close()