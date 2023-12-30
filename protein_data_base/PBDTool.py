import sys
import math

# Global variables
atom_record = []
atom_freq = {}
res_freq = {}
first_model_processed = False

def pdb( path ):
  global first_model_processed
  atom_record = []

  try: 
    with open( path, 'r' ) as file:
      for line in file:
        if line.startswith( 'MODEL' ):
          if first_model_processed:
            break
          else:
            first_model_processed = True

        elif line.startswith( 'ATOM' ):
          atom_serial = int( line[6:11] )
          residue = line[17:20].strip()
          chain = line[21]
          residue_sequence = int( line[22:26] )
          x = float( line[31:38] )
          y = float( line[39:46] )
          z = float( line[47:54] )
          occupancy = float( line[55:60] )
          temp_factor = line[61:66]
          element_name = line[77:78]

          atom_record.append({
              'Atom Serial Number': atom_serial, 
              'Residue Name': residue, 
              'Chain ID': chain, 
              'Residue Sequence': residue_sequence, 
              'x': x,
              'y': y,
              'z': z,
              'Occupany': occupancy, 
              'temp_factor': temp_factor,
              'Element Name': element_name
          })
  
  except FileNotFoundError:
    print( 'The file is not found: {path}.  ' )
    exit(1)

  return atom_record


# Note: Calls Global Variable
def run( c ):
    c = c.strip()

    if c == "help":
      str = "\nList of commands:"
      str += "\n1. help"
      str += "\n2. atomfreq"
      str += "\n3. resfreq"
      str += "\n4. reslength <res_name> <chain_id> <res_seq_num>"
      str += "\n5. tempcheck <decimal: 0.00 - 100.00>"
      str += "\n6. occupancy <decimal: 0.00 - 1.00>"
      str += "\n7. quit"
      print(str)

    elif c == 'atomfreq':
      global atom_freq
      for element in sorted( atom_freq.keys() ):
        print( f' {element}: {atom_freq[element]} ' )

    elif c == 'resfreq':
      global res_freq
      for residue in sorted(res_freq.keys()):
        print( f' {residue}: {res_freq[residue]} ' )

    elif c.startswith( 'reslength' ):
      args = c.split(' ')
      
      # argument validation
      if len( args ) != 4:
        if len( args ) == 0:
          print( 'Missing arguments to reslength. Please use "help" command to learn more detailed about using reslength. ')
        else:
          print( 'Incorrect number of arguments to reslength. ' )
        return 

      res_name = args[1]
      chain_id = args[2]
      res_seq_num = args[3]

      # argument validation
      if len( res_name ) != 3 or not res_name.isupper():
        print("Usage: reslength <res_name> <chain_id> <res_seq_num>. For details about the reslength command, use the 'help' command.")
        return
      if len( chain_id ) != 1 or not chain_id.isupper():
        print("Usage: reslength <res_name> <chain_id> <res_seq_num>. For details about the reslength command, use the 'help' command.")
        return
      if not res_seq_num.isdigit():
        print("Usage: reslength <res_name> <chain_id> <res_seq_num>. For details about the reslength command, use the 'help' command.")
        return

      list_of_res = []
      global atom_record
      for atom in atom_record:
        if atom['Residue Sequence'] == int(res_seq_num) and atom['Residue Name'] == res_name and atom['Chain ID'] == chain_id:
          list_of_res.append(atom)
      
      if len(list_of_res) == 0:
        print("No residue present.")

      distance = 0
      for i in range( len(list_of_res) ):
        for j in range( i + 1, len(list_of_res) ):
          d = math.sqrt( (list_of_res[i]['x'] - list_of_res[j]['x'])**2 + (list_of_res[i]['y'] - list_of_res[j]['y'])**2 + (list_of_res[i]['z'] - list_of_res[j]['z'])**2 )
          if d > distance:
            distance = d

      print( f"{res_name} with sequence number {res_seq_num} in chain {chain_id} has {distance:.2f} angstroms." )
      

    elif c.startswith("tempcheck"):
      args = c.split(" ")
      try:
          threshold = float(args[1])
          if not (0.00 <= threshold <= 100.00):
              raise ValueError
      
      except ValueError:
            print("Usage: tempcheck <decimal>. \n For details about the tempcheck command, use the 'help' command.")     
            return

      for temp in atom_record:
        below_count = sum( 1 for atom in atom_record if float(atom['temp_factor']) < threshold )
        at_count = sum( 1 for atom in atom_record if float(atom['temp_factor']) == threshold )
        above_count = sum( 1 for atom in atom_record if float(atom['temp_factor']) > threshold ) 

      total_atoms = len(atom_record)

      print( f"Temperature factor below {threshold:.2f}: {below_count} / {total_atoms} ({(below_count/total_atoms)*100:.1f}%)")
      print( f"Temperature factor at {threshold:.2f}: {at_count} / {total_atoms} ({(at_count/total_atoms)*100:.1f}%)")
      print( f"Temperature factor above {threshold:.2f}: {above_count} / {total_atoms} ({(above_count/total_atoms)*100:.1f}%)")

    
    elif c.startswith("occupancy"):
      args = c.split(" ")

      try:
          threshold = float(args[1])
          if not (0.00 <= threshold <= 1.00):
              raise ValueError
      
      except ValueError:
            print("Usage: occupancy <decimal>. For details about the occupancy command, use the 'help' command.")     
            return

      # Count atoms below, above, and at the specified occupancy value
      below_count = sum(1 for atom in atom_record if float(atom['Occupany']) < threshold)
      at_count = sum(1 for atom in atom_record if float(atom['Occupany']) == threshold)
      above_count = sum(1 for atom in atom_record if float(atom['Occupany']) > threshold)

      total_atoms = len(atom_record)

      print(f"Occupancy below {threshold:.2f}: {below_count} / {total_atoms} ({(below_count/total_atoms)*100:.1f}%)")
      print(f"Occupancy at {threshold:.2f}: {at_count} / {total_atoms} ({(at_count/total_atoms)*100:.1f}%)")
      print(f"Occupancy above {threshold:.2f}: {above_count} / {total_atoms} ({(above_count/total_atoms)*100:.1f}%)")

    elif c == "quit":
      sys.exit()
    
    else:
      print("Invalid command: Type 'help' for the list of valid commands.")

# Returns a list of stat dicts
def stats(record):
  atom_freq = {}
  res_freq = {}

  for atom in record:
    if atom["Element Name"] in atom_freq.keys():
      atom_freq[atom["Element Name"]] = atom_freq[atom["Element Name"]] + 1
    else:
      atom_freq[atom["Element Name"]] = 1

    if atom["Residue Name"] in res_freq.keys():
      res_freq[atom["Residue Name"]] = res_freq[atom["Residue Name"]] + 1
    else:
      res_freq[atom["Residue Name"]] = 1
  return [atom_freq, res_freq]


# ================ Code Run ================ #

"""
$ python3 PBDTool.py 6lu7.pdb

"""

try:

  atom_record = pdb( sys.argv[1] )
  print( "\n Welcome to the pdb program." )
  print( 'To begin, try typing "help" for the list of valid commands. ' )
  print( f"\n{len(atom_record)} atoms recorded." )

  # Populate information
  results = stats(atom_record)
  atom_freq = results[0]
  res_freq = results[1]
  
  while True:
    print("\nEnter a command:", end=" ")
    command = input()
    run(command)

except IndexError:
  print( "You are missing a PDB file. \n Please enter a file to start." )
