# Brian Enwonwu
# A program that determines the proportions of convs) a column graph showing the proportion of identical, conserved, and variable sites in each set
from itertools import combinations
import math, os
import exceptions, xlsxwriter

#Sample protein sequences below
"""
BAC_FIRMI = 'KPHVNIGTIGHVDHGKTTLTAAISKVLAEK--EGKTATDFAEIDNAPEEKERGITINTSHIEYETDARHYAHIDAPGHADYVKNMITGAAQMDGAILVVAATDGPMPQTREHILLARQVGVEHLVVFLNKTDLVDDEELVDLVEMEVRELLSEYDFPGDDIPVIKGSALKALE----------GDKEQE--QVILDLMKAVDEYIPTPAREDDKPFLMPVEDVFTITGRGTVASGRVDRGVLTTGTEVEIVGL--KDEIQKTTVTGIEMFRKTLDEAQAGDNIGALLRGIDRDQIERGQVLAAPGSIKTHKKFKAEVYVLTKEEGGRHTPFFTNYRPQFYFHTTDVTGVVEL-----PEGVEMVMPGDQVTFDVELIAPVAIEKGLKFTVREGGHTVGAG'
BAC_ENTERO = 'KPHVNVGTIGHVDHGKTTLTAAITTVLAKT--YGGSARAFDQIDNAPEEKARGITINTSHVEYDTPTRHYAHVDCPGHADYVKNMITGAAQMDGAILVVAATDGPMPQTREHILLGRQVGVPFIIVFMNKCDMVDDEELLELVEMEVRELLSAYDFPGDDLPVVRGSALKALE----------GEAEWE--AKIIELAGYLDSYIPEPERAIDKPFLLPIEDVFSISGRGTVVTGRVERGIVKVGEEVEIVGI--K-DTVKSTCTGVEMFRKLLDEGRAGENVGVLLRGIKREDIERGQVLAKPGSIKPHTQFDSEVYILSKDEGGRHTPFFKGYRPQFYFRTTDVTGTIEL-----PEGVEMVMPGDNVNMKVTLIHPIAMDDGLRFAIREGGRTVGAG'
BAC_GCPOS = 'KPHVNIGTIGHVDHGKTTLTAAITKVLHDKFPDLNETKAFDQIDNAPEERQRGITINIAHVEYQTDKRHYAHVDAPGHADYIKNMITGAAQMDGAILVVAATDGPMPQTREHVLLARQVGVPYILVALNKADAVDDEELLELVEMEVRELLAAQEFD-EDAPVVRVSALKALE----------GDAKWV--ASVEELMNAVDESIPDPVRETDKPFLMPVEDVFTITGRGTVVTGRVERGVINVNEEVEIVGI--RPSTTKTTVTGVEMFRKLLDQGQAGDNVGLLLRGVKREDVERGQVVTKPGTTTPHTEFEGQVYILSKDEGGRHTPFFNNYRPQFYFRTTDVTGVVTL-----PEGTEMVMPGDNTNISVKLIQPVAMDEGLRFAIREGGRTVGAG'
BAC_CYANO = 'KPHINIGTVGHVDHGKTTLTAAITMTLAAM--GQAVAKGYDQIDNAPEEKARGITINTAHVEYETANRHYAHVDCPGHADYVKNMITGAAQMDGGILVVAATDGPMPQTREHILLAKQVGVPSLVVFLNKEDMVDDEELLELVELEVRELLSSYEFDGDNIPVIKGSGLQALEAMTKNPKTQRGENPWV--DKIYELMDAVDSYIPTPERDVDKPFLMAVEDVFSITGRGTVATGRIERGVVKVGDNVELVGI--K-DTRATTVTGIEMFKKSLDQGMAGDNAGVLLRGIQKADIERGMVIAKPGSITPHTQFEGEVYVLTDKEGGRKTPFFAGYRPQFYVRTTDVTGTIKAFTSDDGKDVEMVMPGDRIKVTVELINAIAIEQGMRFAIREGGRTIGAG'
"""

valid_amino_acid = ['I', 'V', 'L', 'F', 'C', 'M', 'A', 'G', 'T', 'S', 'W', 'Y', 'P', 'H', 'E', 'Q', 'D', 'N', 'K', 'R']
sequenceList = [] #Stores user-inputted protein sequences

print("Input protein sequences one by one and press Enter key. ")
print("When finished type Done and press Enter key")
i = 1
while True:
  try:
    input_seq = input("Sequence " + str(i) + ": ")
    if input_seq == "Done":
      break
    for char in input_seq:
      if char not in valid_amino_acid and char != '-':
        raise exceptions.Invalid_Character_Exception       
  except exceptions.Invalid_Character_Exception:
    print("Invalid Character found in input. Retry")
    continue
  except exceptions.Empty_Input_Exception:
    print("Input cannot be blank. Retry")
    continue
  sequenceList.append(input_seq)  
  i += 1  

truncated_length = len(sequenceList[0])
for sequence in sequenceList:
  if len(sequence) < truncated_length:
    truncated_length = len(sequence)
    
"""
Function that calculates the probabilties of observed amino acid pairings in the list of protein sequences
@param - sequenceList the list of protein sequences
@return - observed_frequency_dictionary: A dictionary with observed amino acid pairs as keys and their observed probabilities as values
"""    
def compute_observed_probability(sequenceList):
  total_combo_counter = 0
  observed_dictionary = {}
  for index in range(truncated_length):
    temp_list = [] 
    for sequence in sequenceList:
      if sequence[index] != '-':
        temp_list.append(sequence[index])
    temp_list = sorted(temp_list)
    combo_list = list(combinations(temp_list,2))
    total_combo_counter += len(combo_list)

    for i in combo_list:
      observed_dictionary[i] = observed_dictionary.get(i, 0) + 1
  observed_frequency_dictionary = {k: v / total for total in (sum(observed_dictionary.values()),) for k, v in observed_dictionary.items()}
  return observed_frequency_dictionary
    
"""
Function that calculates the probabilities of individual amino acids occurring in the list of protein sequences
@param - sequenceList the list of protein sequences
@return - individual_frequency_dictionary: A dictionary with amino acids as keys and their probabilities as values
"""         
def compute_individual_ultimate(sequenceList):
  individual_frequency_dictionary = {}
  for index in range(truncated_length):
    for sequence in sequenceList:
      if sequence[index] != '-':
        individual_frequency_dictionary[sequence[index]] = individual_frequency_dictionary.get(sequence[index], 0) + 1
  individual_frequency_dictionary = {k: v / total for total in (sum(individual_frequency_dictionary.values()),) for k, v in individual_frequency_dictionary.items()}      
  return individual_frequency_dictionary  

"""
Function that calculates the probabilities of possible amino acid pairings in the list of protein sequences, based on the individual amino acid probabilities
@param - sequenceList the list of protein sequences
@return - expected_frequency_dictionary: A dictionary with expected amino acids as keys and their probabilities as values
"""   
def compute_expected_probability(arg_dict):
  expected_frequency_dictionary = {}
  temp_list = sorted(valid_amino_acid) 
  combo_list = list(combinations(temp_list,2))
  for amino in valid_amino_acid:
    temp_tuple = (amino, amino)
    combo_list.append(temp_tuple)
  
  #print(combo_list)
  for i in combo_list:
    if i[0] == i[1]:
      expected_frequency_dictionary[i] = arg_dict.get(i[0]) * arg_dict.get(i[1])
    elif i[0] != i[1]:
      expected_frequency_dictionary[i] = 2 * arg_dict.get(i[0]) * arg_dict.get(i[1])
  return expected_frequency_dictionary

"""
Function that calculates log odds probabilities of amino acid pairs from the observed ad expected probabilities
@param - obvs_dict the list of protein sequences
@param - expected_dict: the list of protein sequences
@return - log_dictionary: A dictionary with expected amino acids as keys and their log odds probabilities as values
"""     
def compute_log(obvs_dict, expected_dict):
  log_dictionary = {}
  for k,v in expected_dict.items():
    if k in obvs_dict:
      log_dictionary[k] = round(math.log2(obvs_dict.get(k)/ v ) * 2)
    else:
      log_dictionary[k] = 0
  return log_dictionary    

"""
Function that finds the two most similar sequences in a list of protein sequences
@param - sequence_list: the list of protein sequences
@return - similar_list: a list containing the two most similar sequences 
"""
def find_most_similar(sequence_list):  
  sorted_list = sequence_list
  combo_list = list(combinations(sorted_list,2))
  highest_similarity_score = 0
  sim_string_one = ""
  sim_string_two = ""
  similar_list = []
  for string_pair in combo_list:
    current_similarity_score = 0
    for i in range(len(string_pair[0])):
      for x in range(len(string_pair[0][i])):
        for y in range(len(string_pair[1][i])):
          if string_pair[0][i][x] == string_pair[1][i][y]:
            current_similarity_score += 1
       
    if current_similarity_score > highest_similarity_score:
      highest_similarity_score = current_similarity_score 
      sim_string_one =  string_pair[0]  
      sim_string_two =  string_pair[1] 
  similar_list.append(sim_string_one)
  similar_list.append(sim_string_two)     
  return similar_list

GAP_PENALTY = -3 # Used in alignment matrix

"""
Function that creates an alignment matrix from two protein sequences
@param - sequence_one: the first protein sequence
@param - sequence_one: the second protein sequence
@param - log_odds_dictionary: dictionary containing amino acid pairs and their log odds probabilities
@return - similar_list: a list containing the two most similar sequences 
"""
def create_alignment_matrix(sequence_one, sequence_two, log_odds_dictionary):
  # Initialize the matrix with zeros
  matrix = [[0] * (len(sequence_two) + 1) for _ in range(len(sequence_one) + 1)]

  # Initialize the first column and row with gap penalties
  for i in range(1, len(sequence_one) + 1):
    matrix[i][0] = matrix[i - 1][0] + GAP_PENALTY

  for j in range(1, len(sequence_two) + 1):
    matrix[0][j] = matrix[0][j - 1] + GAP_PENALTY

  # Fill the matrix based on match/mismatch scores and gap penalties, DOESNT WORK THO FIX THIS WHEN U AWAKE FAM
  for i in range(1, len(sequence_one) + 1):
    for j in range(1, len(sequence_two) + 1):
      list_of_max_values = []
      for x in range(len(sequence_one[i - 1])):
        for y in range(len(sequence_two[j - 1])):
          temp_tuple = (sequence_one[i - 1][x], sequence_two[j - 1][y])
          diagonal_value = matrix[i - 1][j - 1] + log_odds_dictionary.get(temp_tuple, 0) #initialize to max_diagonal_value then replace
          list_of_max_values.append(diagonal_value)
             
          coming_from_above_value = matrix[i - 1][j] + GAP_PENALTY
          list_of_max_values.append(coming_from_above_value)
             
          coming_from_left_value = matrix[i][j - 1] + GAP_PENALTY
          list_of_max_values.append(coming_from_left_value)
                               
          matrix[i][j] = max(list_of_max_values)
  return matrix  

indi = compute_individual_ultimate(sequenceList)
obvs = compute_observed_probability(sequenceList)
expec = compute_expected_probability(indi)
log_odds_values = compute_log(obvs, expec)
  
sequence_test_list = find_most_similar(sequenceList)  #finds two of the most aligned
sequence_one = sequence_test_list[0]
sequence_two = sequence_test_list[1]

alignment_matrix = create_alignment_matrix(sequence_one,sequence_two, log_odds_values) #ad dsomething to print out alignment matrix

#Writes alignment matrix to excel sheet
script_dir = os.path.dirname(os.path.abspath(__file__))
file_path = os.path.join(script_dir, 'matrix.xlsx')
with xlsxwriter.Workbook(file_path) as workbook:
    worksheet = workbook.add_worksheet()
    for row_num, data in enumerate(alignment_matrix):
      worksheet.write_row(row_num, 0, data)
