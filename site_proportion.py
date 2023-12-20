# Brian Enwonwu
# A program that determines the proportions of  identical, conserved, and variable sites in a set of protein sequences based on Grantham Distances and plots a pie chart of the data
from itertools import *
import numpy as np
import matplotlib.pyplot as plt
import exceptions

"""
Samples to test program
ANIM_ARTH = 'QAHCNVGTIGHVDHGKTTLTAAITKIQSGK--GLADYVSYDQIDRAPEEKARGITINACHIGYATKERTYAHTDCPGHADYIKNMISGASQMDGAILVVAATDGQMPQTREHLLLAKQVGIERIVVFINKADLV-DQEVLELVEIEMREMLTDFGFDGVNSPVICGSALLALR----------GDQSSFGVPAIEQLLQHCDNYIPTPKRDTTAPFILPIDNAFTVPGRGTVVVGTIKRGTILRNADADLLGF---NQNLKTTVSDIQIFRKSVPQALAGENVGALLRGIKISAVERGMLLCATGSENISNHFEASMYLLSRAEGGRVKPMLSKYIQQLFSMTWNTPARIDM-----VPHESMLMPGEHSKVRVTLMRKMVMTPGQAFTIRENGMTVATG'
ANIM_CHORD = 'KPHVNVGTIGHVDHGKTTLTAAITKILAEG--GGAKFKKYEEIDNAPEERARGITINAAHVEYSTAARHYAHTDCPGHADYVKNMITGTAPLDGCILVVAANDGPMPQTREHLLLARQIGVEHVVVYVNKADAVQDSEMVELVELEIRELLTEFGYKGEETPVIVGSALCALE----------GRDPELGLKSVQKLLDAVDTYIPVPARDLEKPFLLPVEAVYSVPGRGTVVTGTLERGILKKGDECELLGH---SKNIRTVVTGIEMFHKSLERAEAGDNLGALVRGLKREDLRRGLVMVKPGSIKPHQKVEAQVYILSKEEGGRHKPFVSHFMPVMFSLTWNMACRIIL-----PPEKELAMPGEDLKFNLILRQPMILEKGQRFTLRDGNRTIGTG'
ANIM_CNID = 'KPHINIGTIGHVDHGKTTLTAAITKVLSEK--GGSKFKDYADIDNAPEERARGITINASHVEYETDTRHYGHIDCPGHADYIKNMITGAAQMDGAILVVAATDGQMPQTREHLLLANQIGVKNLCVFINKADMVDDKEIMDLVEMEIRELLTEYGYDGDNTPVIGGSALCALE----------GKKPELGVQKIQELLAAVDSHIPLPKRDLDKPFLMPVEDSFSISGRGTVITGSIERGIVKKGDELELVGH--SNVPIKTVATGLEMFHKSLEQGQAGDNLGALVRGLKREDVKRGMVLCAPGTVKAYTKCKAQVYILKKEEGGRHKPFVSNYTPQMYVRTGDVAATITL-----DAGKEFVMPGEDASFSLTLMHPTPLEKGLRFTMREGSKTVGTG'
ANIM_SPONG = 'KPHVNVGTIGHVDHGKTTLTAAITRVLAEE--GGAAFRGMDEIDNAPEEKARGVTIAIAHVEYETKERHYAHVDCPGHADYIKNMITGAAQMDGAVLVVSAPDGPMPQTREHILLARQVEVPSMVVALNKVDMMDDEELLELVELEIRELLSRYEFPGDDTPIVRVSGLKALD----------GDPEAM--QGVRELVQTMDDYIPLPTRITDQAFLMPIEDVFGIKGRGTVVTGRVERGQLSVGQEVEIIRS--G-DVRKTVATGLEMFHKLLDTTEAGDAVGVLLRGVDRDEVERGQVLVAPGSMRPYRQAEAEVYVLSQQEGGRHTPFFTGYKPQFYIRTSDITGEIGL-----PEGVEMVMPGDNITMTVNLITPIAIEEGLRFAIREGGRTVGAG'
"""
valid_amino_acid = ['I', 'V', 'L', 'F', 'C', 'M', 'A', 'G', 'T', 'S', 'W', 'Y', 'P', 'H', 'E', 'Q', 'D', 'N', 'K', 'R']
sequenceList = []

print("Input protein sequences one by one and press Enter key. ")
print("When finished type Complete and press Enter key")
i = 1
while True:
  try:
    input_seq = input("Sequence " + str(i) + ": ")
    if input_seq == "Complete":
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

grantham_dictionary = { ('A', 'A'): '0',  ('A', 'C'): '195',  ('A', 'D'): '126', ('A', 'E'): '107',  ('A', 'F'): '113',  ('A', 'G'): '60',  ('A', 'H'): '86',  ('A', 'I'): '94',  ('A', 'K'): '106',  ('A', 'L'): '96',  ('A', 'M'): '84',
  ('A', 'N'): '111',  ('A', 'P'): '27',  ('A', 'Q'): '91',  ('A', 'R'): '112',  ('A', 'S'): '99',  ('A', 'T'): '58',  ('A', 'V'): '64',  ('A', 'W'): '148',  ('A', 'Y'): '112',  ('C', 'A'): '195',  ('C', 'C'): '0',  ('C', 'D'): '154',
  ('C', 'E'): '170',  ('C', 'F'): '205',  ('C', 'G'): '159',  ('C', 'H'): '174',  ('C', 'I'): '198',  ('C', 'K'): '202',  ('C', 'L'): '198',  ('C', 'M'): '196',  ('C', 'N'): '139',  ('C', 'P'): '169',  ('C', 'Q'): '154',
  ('C', 'R'): '180',  ('C', 'S'): '112',  ('C', 'T'): '149',  ('C', 'V'): '192',  ('C', 'W'): '215',  ('C', 'Y'): '194',  ('D', 'A'): '126',  ('D', 'C'): '154',  ('D', 'D'): '0',  ('D', 'E'): '45',  ('D', 'F'): '177',  
  ('D', 'G'): '94',  ('D', 'H'): '81',  ('D', 'I'): '168',  ('D', 'K'): '101',  ('D', 'L'): '172',  ('D', 'M'): '160',  ('D', 'N'): '23',  ('D', 'P'): '108',  ('D', 'Q'): '61',  ('D', 'R'): '96',  ('D', 'S'): '65',  
  ('D', 'T'): '85',  ('D', 'V'): '152',  ('D', 'W'): '181',  ('D', 'Y'): '160',  ('E', 'A'): '107',  ('E', 'C'): '170',  ('E', 'D'): '45',  ('E', 'E'): '0',  ('E', 'F'): '140',  ('E', 'G'): '98',  ('E', 'H'): '40',  
  ('E', 'I'): '134',  ('E', 'K'): '56',  ('E', 'L'): '138',  ('E', 'M'): '126',  ('E', 'N'): '42',  ('E', 'P'): '93',  ('E', 'Q'): '29',  ('E', 'R'): '54',  ('E', 'S'): '80',  ('E', 'T'): '65',  ('E', 'V'): '121',  
  ('E', 'W'): '152',  ('E', 'Y'): '122',  ('F', 'A'): '113',  ('F', 'C'): '205',  ('F', 'D'): '177',  ('F', 'E'): '140',  ('F', 'F'): '0',  ('F', 'G'): '153',  ('F', 'H'): '100',  ('F', 'I'): '21',  ('F', 'K'): '102',
  ('F', 'L'): '22',  ('F', 'M'): '28',  ('F', 'N'): '158',  ('F', 'P'): '114',  ('F', 'Q'): '116',  ('F', 'R'): '97',  ('F', 'S'): '155',  ('F', 'T'): '103',  ('F', 'V'): '50',  ('F', 'W'): '40',  ('F', 'Y'): '22',  
  ('G', 'A'): '60',  ('G', 'C'): '159',  ('G', 'D'): '94',  ('G', 'E'): '98',  ('G', 'F'): '153',  ('G', 'G'): '0',  ('G', 'H'): '98',  ('G', 'I'): '135',  ('G', 'K'): '127',  ('G', 'L'): '138',  ('G', 'M'): '127',  
  ('G', 'N'): '80',  ('G', 'P'): '42',  ('G', 'Q'): '87',  ('G', 'R'): '125',  ('G', 'S'): '56',  ('G', 'T'): '59',  ('G', 'V'): '109',  ('G', 'W'): '184',  ('G', 'Y'): '147',
  ('H', 'A'): '86',  ('H', 'C'): '174',  ('H', 'D'): '81',  ('H', 'E'): '40',  ('H', 'F'): '100',  ('H', 'G'): '98',
  ('H', 'H'): '0',  ('H', 'I'): '94',  ('H', 'K'): '32',  ('H', 'L'): '99',  ('H', 'M'): '87',  ('H', 'N'): '68',  ('H', 'P'): '77',  ('H', 'Q'): '24',
  ('H', 'R'): '29',  ('H', 'S'): '89',  ('H', 'T'): '47',  ('H', 'V'): '84',  ('H', 'W'): '115',  ('H', 'Y'): '83',
  ('I', 'A'): '94',  ('I', 'C'): '198',  ('I', 'D'): '168',  ('I', 'E'): '134',  ('I', 'F'): '21',  ('I', 'G'): '135',
  ('I', 'H'): '94',  ('I', 'I'): '0',  ('I', 'K'): '102',  ('I', 'L'): '5',  ('I', 'M'): '10',  ('I', 'N'): '149',  ('I', 'P'): '95',  ('I', 'Q'): '109',
  ('I', 'R'): '97',  ('I', 'S'): '142',  ('I', 'T'): '89',  ('I', 'V'): '29',  ('I', 'W'): '61',  ('I', 'Y'): '33',  ('K', 'A'): '106',  ('K', 'C'): '202',
  ('K', 'D'): '101',  ('K', 'E'): '56',  ('K', 'F'): '102',  ('K', 'G'): '127',  ('K', 'H'): '32',  ('K', 'I'): '102',
  ('K', 'K'): '0',  ('K', 'L'): '107',  ('K', 'M'): '95',  ('K', 'N'): '94',  ('K', 'P'): '103',  ('K', 'Q'): '53',
  ('K', 'R'): '26',  ('K', 'S'): '121',  ('K', 'T'): '78',  ('K', 'V'): '97',  ('K', 'W'): '110',  ('K', 'Y'): '85',  ('L', 'A'): '96',  ('L', 'C'): '198',  ('L', 'D'): '172',  ('L', 'E'): '138',  ('L', 'F'): '22',
  ('L', 'G'): '138',  ('L', 'H'): '99',  ('L', 'I'): '5',  ('L', 'K'): '107',  ('L', 'L'): '0',
  ('L', 'M'): '15',  ('L', 'N'): '153',  ('L', 'P'): '98',  ('L', 'Q'): '113',  ('L', 'R'): '102',  ('L', 'S'): '145',  ('L', 'T'): '92',  ('L', 'V'): '32',  ('L', 'W'): '61',
  ('L', 'Y'): '36',  ('M', 'A'): '84',  ('M', 'C'): '196',  ('M', 'D'): '160',  ('M', 'E'): '126',  ('M', 'F'): '28',  ('M', 'G'): '127',  ('M', 'H'): '87',
  ('M', 'I'): '10',  ('M', 'K'): '95',  ('M', 'L'): '15',  ('M', 'M'): '0',  ('M', 'N'): '142',  ('M', 'P'): '87',
  ('M', 'Q'): '101',  ('M', 'R'): '91',  ('M', 'S'): '135',  ('M', 'T'): '81',  ('M', 'V'): '21',  ('M', 'W'): '67',  ('M', 'Y'): '36',
  ('N', 'A'): '111',  ('N', 'C'): '139',  ('N', 'D'): '23',  ('N', 'E'): '42',  ('N', 'F'): '158',  ('N', 'G'): '80',  ('N', 'H'): '68',  ('N', 'I'): '149',  ('N', 'K'): '94',  ('N', 'L'): '153',  ('N', 'M'): '142',  ('N', 'N'): '0',  ('N', 'P'): '91',  ('N', 'Q'): '46',  ('N', 'R'): '86',
  ('N', 'S'): '46',  ('N', 'T'): '65',  ('N', 'V'): '133',  ('N', 'W'): '174',  ('N', 'Y'): '143',  ('P', 'A'): '27',  ('P', 'C'): '169',  ('P', 'D'): '108',  ('P', 'E'): '93',  ('P', 'F'): '114',  ('P', 'G'): '42',
  ('P', 'H'): '77',  ('P', 'I'): '95',  ('P', 'K'): '103',  ('P', 'L'): '98',  ('P', 'M'): '87',  ('P', 'N'): '91',  ('P', 'P'): '0',  ('P', 'Q'): '76',  ('P', 'R'): '103',  ('P', 'S'): '74',
  ('P', 'T'): '38',  ('P', 'V'): '68',  ('P', 'W'): '147',  ('P', 'Y'): '110',  ('Q', 'A'): '91',  ('Q', 'C'): '154',  ('Q', 'D'): '61',
  ('Q', 'E'): '29',  ('Q', 'F'): '116',  ('Q', 'G'): '87',  ('Q', 'H'): '24',  ('Q', 'I'): '109',  ('Q', 'K'): '53',
  ('Q', 'L'): '113',  ('Q', 'M'): '101',  ('Q', 'N'): '46',  ('Q', 'P'): '76',  ('Q', 'Q'): '0',  ('Q', 'R'): '43',  ('Q', 'S'): '68',  ('Q', 'T'): '42',  ('Q', 'V'): '96',
  ('Q', 'W'): '130',  ('Q', 'Y'): '99',  ('R', 'A'): '112',  ('R', 'C'): '180',  ('R', 'D'): '96',  ('R', 'E'): '54',  ('R', 'F'): '97',  ('R', 'G'): '125',
  ('R', 'H'): '29',  ('R', 'I'): '97',  ('R', 'K'): '26',  ('R', 'L'): '102',  ('R', 'M'): '91',  ('R', 'N'): '86',  ('R', 'P'): '103',  ('R', 'Q'): '43',  ('R', 'R'): '0',  ('R', 'S'): '110',  ('R', 'T'): '71',
  ('R', 'V'): '96',  ('R', 'W'): '101',  ('R', 'Y'): '77',  ('S', 'A'): '99',  ('S', 'C'): '112',  ('S', 'D'): '65',  ('S', 'E'): '80',  ('S', 'F'): '155',
  ('S', 'G'): '56',  ('S', 'H'): '89',  ('S', 'I'): '142',  ('S', 'K'): '121',  ('S', 'L'): '145',  ('S', 'M'): '135',  ('S', 'N'): '46',  ('S', 'P'): '74',  ('S', 'Q'): '68',
  ('S', 'R'): '110',  ('S', 'S'): '0',  ('S', 'T'): '58',  ('S', 'V'): '124',  ('S', 'W'): '177',  ('S', 'Y'): '144',  ('T', 'A'): '58',  ('T', 'C'): '149',  ('T', 'D'): '85',  ('T', 'E'): '65',  ('T', 'F'): '103',  ('T', 'G'): '59',  ('T', 'H'): '47',  ('T', 'I'): '89',  ('T', 'K'): '78',
  ('T', 'L'): '92',  ('T', 'M'): '81',  ('T', 'N'): '65',  ('T', 'P'): '38',  ('T', 'Q'): '42',  ('T', 'R'): '71',  ('T', 'S'): '58',  ('T', 'T'): '0',
  ('T', 'V'): '69',  ('T', 'W'): '128',  ('T', 'Y'): '92',  ('V', 'A'): '64',  ('V', 'C'): '192',  ('V', 'D'): '152',  ('V', 'E'): '121',
  ('V', 'F'): '50',  ('V', 'G'): '109',  ('V', 'H'): '84',  ('V', 'I'): '29',  ('V', 'K'): '97',  ('V', 'L'): '32',  ('V', 'M'): '21',
  ('V', 'N'): '133',  ('V', 'P'): '68',  ('V', 'Q'): '96',  ('V', 'R'): '96',  ('V', 'S'): '124',  ('V', 'T'): '69',  ('V', 'V'): '0',  ('V', 'W'): '88',  ('V', 'Y'): '55',
  ('W', 'A'): '148',  ('W', 'C'): '215',  ('W', 'D'): '181',  ('W', 'E'): '152',  ('W', 'F'): '40',  ('W', 'G'): '184',  ('W', 'H'): '115',
  ('W', 'I'): '61',  ('W', 'K'): '110',  ('W', 'L'): '61',  ('W', 'M'): '67',  ('W', 'N'): '174',  ('W', 'P'): '147',  ('W', 'Q'): '130',
  ('W', 'R'): '101',  ('W', 'S'): '177',  ('W', 'T'): '128',  ('W', 'V'): '88',  ('W', 'W'): '0',  ('W', 'Y'): '37',  ('Y', 'A'): '112',  ('Y', 'C'): '194',
  ('Y', 'D'): '160',  ('Y', 'E'): '122',  ('Y', 'F'): '22',  ('Y', 'G'): '147',  ('Y', 'H'): '83', ('Y', 'I'): '33',  ('Y', 'K'): '85',  ('Y', 'L'): '36',
  ('Y', 'M'): '36',  ('Y', 'N'): '143',  ('Y', 'P'): '110',  ('Y', 'Q'): '99',  ('Y', 'R'): '77',  ('Y', 'S'): '144',  ('Y', 'T'): '92',  ('Y', 'V'): '55',
  ('Y', 'W'): '37',  ('Y', 'Y'): '0'}

conserved_count = 0
variable_count = 0
identical_count = 0
for index in range(truncated_length):
  amino_acid_list = []
  max_score = 0
  for sequence in sequenceList:
    if sequence[index] != '-':
      amino_acid_list.append(sequence[index])
    combination_list = list(combinations(amino_acid_list, 2)) 
    for element in combination_list:
      current_score = int(grantham_dictionary.get(element, 0)) 
      if current_score > max_score:
        max_score = current_score
  if max_score == 0:
    identical_count += 1
  elif max_score >= 100:
    variable_count += 1
  else:
    conserved_count += 1     
        
print("Conserved count of Phylo is "+ str(100 * conserved_count / truncated_length))    
print("Variable count of Phylo is "+ str(100 * variable_count / truncated_length))  
print("Indentical count of Phylo is "+ str(100 * identical_count/ truncated_length)) 

#Pie chart creation
pie_chart = np.array([conserved_count, variable_count, identical_count])
mylabels = ["Conserved", "Variable", "Identical"]
plt.pie(pie_chart, labels = mylabels, autopct='%1.1f%%')
plt.show()