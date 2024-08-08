import copy

from flask import Flask, render_template, request
import itertools
# Configure application
app = Flask(__name__)

# Ensure templates are auto-reloaded
app.config["TEMPLATES_AUTO_RELOAD"] = True
# -*- coding: utf-8 -*-

def hl_explain(c1, c2, color, explanation):
  # If a function impacts the possibilities array at cell (c1, c2), it calls this function to document the change.
  global something_changed
  something_changed = True
  # Add a length-2 array at cell (c1, c2) which records the explanation and the color
  explanations[c1][c2].append([explanation,color])

# Possibilities is initialized as an array with 1-9 for each corresponding cell on the sudoku board
# The functions below attempt to narrow down the possibilities array in each cell
poss = [[[1, 2, 3, 4, 5, 6, 7, 8, 9] for j in range(9)] for i in range(9)]
sudoku_board = []
tuples = ["single", "pair", "triple", "quadruple"]
explanations = [[[] for j in range(9)] for i in range(9)]

def naked_single(poss):
  """
  For each cell with a determined value, eliminate that value from
  the possibilities of cells that share a row, column, or block with it
  """
  for n in range(9):
      for m in range(9):
          # run the following for all cells x who have 1 possibility (are already determined)
          if int(len(poss[n][m])) == 1:
              x = poss[n][m][0]
              # go through x's row, and remove x from the possibilities of all the other cells...
              for m2 in range(9):
                  # ... while making sure those other cells are not x itself, and have the possibility in the first place.
                  if m2 != m and x in poss[n][m2]:
                      poss[n][m2].remove(x)
                      hl_explain(n,m2,"r",f"Removed {x} from r{n+1}c{m2+1} possibilities | Naked single in row: r{n+1}c{m+1}")
                      
              # go through x's column, and remove x from the possibilities of all the other cells...
              for n2 in range(9):
                  # ... while making sure those other cells are not x itself, and have the possibility in the first place.
                  if n2 != n and x in poss[n2][m]:
                      poss[n2][m].remove(x)
                      hl_explain(n2,m,"r",f"Removed {x} from r{n2+1}c{m+1} possibilities | Naked single in column: r{n+1}c{m+1}")
                      
              # go through x's block, and remove x from the possibilities of all the other cells...
              for n2 in range(n - n % 3, n - n % 3 + 3):
                  for m2 in range(m - m % 3, m - m % 3 + 3):
                      # ... while making sure those other cells are not x itself, and have the possibility in the first place
                      if (n2 != n or m2 != m) and x in poss[n2][m2]:
                          poss[n2][m2].remove(x)
                          hl_explain(n2,m2,"r",f"Removed {x} from r{n2+1}c{m2+1} possibilities | Naked single in block: r{n+1}c{m+1}")
                          
  return poss


def naked_pair_helper(n, m, n2, m2):
    """ 
    Given the coords of a naked pair, eliminates the values in the pair from the other cells
    in the row/column/block
    """
    # check if the pair share the row
    if n == n2:
        #run through all coordinates in the row that are not the coordinates of the pair
        for m3 in range(9):
            if m3 != m and m3 != m2:
                # run through each value of the pair, checking if the cell has that value and, if so, removing it
                for c in range(2):
                    if poss[n][m][c] in poss[n][m3]:
                        poss[n][m3].remove(poss[n][m][c])
                        hl_explain(n,m3,"r",f"Removed {poss[n][m][c]} from r{n+1}c{m3+1} possibilities | Naked pair in row: r{n+1}c{m+1} and r{n2+1}c{m2+1}")     

    # check if the pair share the column
    if m == m2:
      #run through all coordinates in the column that are not the coordinates of the pair
        for n3 in range(9):
            if n3 != n and n3 != n2:
                # run through each value of the pair, checking if the cell has that value and, if so, removing it
                for c in range(2):
                    if poss[n][m][c] in poss[n3][m]:
                        poss[n3][m].remove(poss[n][m][c])
                        hl_explain(n3,m,"r",f"Removed {poss[n][m][c]} from r{n3+1}c{m+1} possibilities | Naked pair in column: r{n+1}c{m+1} and r{n2+1}c{m2+1}")     

    # check if the pair share the block
    if (n - n % 3 == n2 - n2 % 3) and (m - m % 3 == m2 - m2 % 3):
        #run through all coordinates in the block that are not the coordinates of the pair
        for n3 in range(n - n % 3, n - n % 3 + 3):
            for m3 in range(m - m % 3, m - m % 3 + 3):
                if [m3, n3] != [m2, n2] and [m3, n3] != [m, n]:
                    # run through each value of the pair, checking if the cell has that value and, if so, removing it
                    for c in range(2):
                        if poss[n][m][c] in poss[n3][m3]:
                            poss[n3][m3].remove(poss[n][m][c])
                            hl_explain(n3,m3,"r",f"Removed {poss[n][m][c]} from r{n3+1}c{m3+1} possibilities | Naked pair in block: r{n+1}c{m+1} and r{n2+1}c{m2+1}")
    return poss


def naked_pair(poss):
    """ 
    Finds the coords of naked pair (a pair of cells with exactly two identical possibilities)
    then calls naked_pair_helper to eliminate values from poss
    """
    for n in range(9):
        for m in range(9):
            #run the following for all cells which have two possibilities:
            if int(len(poss[n][m])) == 2:
                # check if there is another cell with the same two possibilities in row
                for m2 in range(9):
                    if m2 != m and poss[n][m2] == poss[n][m]:
                        poss = naked_pair_helper(n,m,n,m2)
                # check if there is another pair with the same two possibilities in column
                for n2 in range(9):
                    if n2 != n and poss[n2][m] == poss[n][m]:
                        poss = naked_pair_helper(n,m,n2,m)
                # check if there is another pair with the same two possibilities in block
                for n2 in range(n - n % 3, n - n % 3 + 3):
                    for m2 in range(m - m % 3, m - m % 3 + 3):
                        if [n2,m2] != [n,m] and poss[n2][m2] == poss[n][m]:
                            poss = naked_pair_helper(n,m,n2,m2)
    return poss

def row_pointer(box, num):
  """
  Takes a row in a box and identifies the coordinates of the row and box
  that it points to
  """
  # compute coordinates of pointing row
  row = [(box[0] + num, box[1] + j) for j in range(0, 3)]
  
  pointed_box = []

  # iterate through box rows
  for i in range(box[0], box[0] + 3):

    # if box row isn't pointing row, add row to pointed_box
    if i != row[0][0]:
      pointed_box += [(i, box[1] + j) for j in range(0, 3)]
    
    # if box row is pointing row, add rest of row to pointed_row
    else:
      full_row = [(i, j) for j in range(0, 9)]
      pointed_row = [cell for cell in full_row if cell not in row]
  
  # return coordinates of the pointer row and the row and box it points to
  return [row, pointed_row, pointed_box]

def col_pointer(box, num):
  """
  Takes a column in a box and identifies the coordinates of the column and
  box that it points to
  """
  # compute coordinates of pointing column
  col = [(box[0] + i, box[1] + num) for i in range(0, 3)]

  pointed_box = []

  # iterate through rows
  for j in range(box[1], box[1] + 3):

    # if box column isn't pointing column, add column to pointed_box
    if j != col[0][1]:
      pointed_box += [(box[0] + i, j) for i in range(0, 3)]
    
    # if box column is pointing column, add rest of column to pointing_column
    else:
      full_col = [(i, j) for i in range(0, 9)]
      pointed_col = [cell for cell in full_col if cell not in col]
  
  # return coordinates of pointer column and the column and box it points to
  return [col, pointed_col, pointed_box]

def pointing_remove(poss, digit, pointed, coords):
  """
  Removes a digit from a set of coordinates that a pointing tuple points to
  """
  # keep track of coords digit was removed from
  removes = []

  # iterate through digits in pointed set
  for i in range(0, 6):

    # remove digit from pointed set
    if digit in pointed[i]:
      poss[coords[i][0]][coords[i][1]].remove(digit)
      removes.append(coords[i])

  return [poss, removes]

def find_format_ptuple(digit, pointer, coords):
  """
  Determine if tuple is a pair or a triple and if so, returns it
  properly formatted
  """
  # initialize list of coordinates in tuple
  ptuple = []

  # find tuple coordinates
  for i in range(0, 3):
    if digit in pointer[i]:
      ptuple.append(coords[i])

  # return false if tuple is length 1
  if len(ptuple) == 1:
    return False

  ptuple_formatted = []

  # format ptuple in r(row)c(column) format
  for coord in ptuple:
    ptuple_formatted.append('r' + str(coord[0]+1) + 'c' + str(coord[1]+1))
  
  # remove brackets and apostrophes
  str_ptuple = str(ptuple_formatted)[1 : -1]
  str_ptuple = str_ptuple.replace("'", "")

  return str_ptuple

def pointing_explain(digit, str_ptuple, removes):
  """
  Explains reason for digit removal calling hl_explain
  """
  # check if any digit was removed
  if len(removes) != 0:
    
    # iterate through cells digit was removed from and determine explanation
    for remove in removes:
      num_tuple = len(str_ptuple.split(","))
      hl_explain(remove[0], remove[1], "y", f"Removed {digit} from r{remove[0]+1}c{remove[1]+1} possibilities | Pointing {tuples[num_tuple - 1]} at {str_ptuple}")

  return

def pointing(pointer_coords, poss):
  """
  Finds pointing tuples and calls pointing_remove to remove possibilities
  """
  # convert pointer coordinates to possibilities in poss
  pointer = [poss[i][j] for i, j in pointer_coords[0]]
  pointed_rowcol = [poss[i][j] for i, j in pointer_coords[1]]
  pointed_box = [poss[i][j] for i, j in pointer_coords[2]]

  # iterate through digits
  for digit in range(1, 10):

    # make sure that pointing tuple is a pair or a triple
    str_ptuple = find_format_ptuple(digit, pointer, pointer_coords[0])
    if str_ptuple != False:

      # if digit isn't in rest of the row/column, remove from pointed box
      if not any(digit in cell for cell in pointed_rowcol):
          poss, removes = pointing_remove(poss, digit, pointed_box, pointer_coords[2])
          pointing_explain(digit, str_ptuple, removes)
    
      # if digit isn't in the rest of box, remove from pointed row/column
      elif not any(digit in cell for cell in pointed_box):
        poss, removes = pointing_remove(poss, digit, pointed_rowcol, pointer_coords[1])
        pointing_explain(digit, str_ptuple, removes)

  return poss

def run_pointing(poss):
  """
  Runs pointing on each box row and column
  """
  # iterate through boxes
  for box in boxes:

    # iterate through rows and columns of box
    for num in range(0, 3):

      # run pointing
      poss = pointing(row_pointer(box, num), poss)
      poss = pointing(col_pointer(box, num), poss)

  return poss
def coordinates(poss):
  """
  Create a 2D array of each digit's possible coordinates on the board
  """
  coords = [[] for i in range(0, 9)]
  # iterate through cells in poss
  for n in range(0, 9):
    for m in range(0, 9):

      # iterate through sudoku digits
      for digit in range(1, 10):

        # add coordinate to digit's list of possible coordinates if digit is a possibility
        if digit in poss[n][m]:
          coords[digit - 1].append((n, m))

  return coords

# calculate coords array
coords = coordinates(poss)

# generate list of combinations for hidden singles, doubles, triples, quadruples
combos = []
# iterate through combination size
for i in range(1, 5):
  # add all combinations
  combos += list(itertools.combinations(range(1, 10), i))


def row_coordinates(n):
  """
  Generate 2D array of possible coordinates per sudoku number for a given row
  """
  row_coords = [[] for i in range(0, 9)]
  # iterate through sudoku number's coordinates
  for i in range(0, 9):
    for coord in coords[i]:

      # add coordinate to array if row number matches
      if coord[0] == n:
        row_coords[i].append(coord)

  return row_coords


def col_coordinates(m):
  """
  Generate 2D array of possible coordinates per sudoku number for a given column
  """
  col_coords = [[] for i in range(0, 9)]
  # iterate through sudoku number's coordinates
  for i in range(0, 9):
    for coord in coords[i]:

      # add coordinate to array if column number matches
      if coord[1] == m:
        col_coords[i].append(coord)

  return col_coords

# initialize top-left coordinate of each sudoku box
boxes = [[n, m] for n in range(0, 9, 3) for m in range(0, 9, 3)]
def box_coordinates(b):
  """
  Generate 2D array of possible coordinates per sudoku number for a given box
  """
  box_coords = [[] for i in range(0, 9)]
  # iterate through sudoku number's coordinates
  for i in range(0, 9):
    for coord in coords[i]:

      # add coordinate to array if coordinate in box
      if coord[0] in range(boxes[b][0], boxes[b][0] + 3) and coord[1] in range(boxes[b][1], boxes[b][1] + 3):
        box_coords[i].append(coord)

  return box_coords


def hidden(poss, group):
  """
  Find hidden subsets in a given row, column, or box

  For each combination of sudoku numbers, find all possible coordinates
  If the number of cells they can be located is equal to the number of sudoku numbers,
  the combination is a hidden subset
  """
  # iterate through combinations of sudoku numbers
  for combo in combos:

    # initialize empty list to put possible coordinates of sudoku numbers in combination
    combo_coords = []

    # iterate through numbers in combination
    for num in combo:

      # make combo_coords the union of combo_coords and the number's coordinates
      if not set(group[num - 1]).issubset(set(combo_coords)):
        for coord in group[num - 1]:
          combo_coords.append(coord)
    
    # check if number of possible coordinates is equal to number of sudoku numbers
    if len(combo_coords) == len(combo):

      # iterate through possible coordinates
      for combo_coord in combo_coords:

        # remove other numbers not in hidden single/double/triple/quadruple
        removes = []
        for digit in poss[combo_coord[0]][combo_coord[1]]:
          if digit not in combo:
            removes.append(digit)
        for digit in removes:
          poss[combo_coord[0]][combo_coord[1]].remove(digit)

          # reformat combo for explanation
          combo_str = str(combo)
          if len(combo) == 1:
            combo_str = combo_str.replace(",", "")
          print(len(combo))
          hl_explain(combo_coord[0], combo_coord[1], "o", f"Removed {digit} from r{combo_coord[0]+1}c{combo_coord[1]+1} possibilities | Hidden {tuples[len(combo)-1]} of {combo_str}")

  return poss
  

def run_hidden(poss):
  """
  Run hidden on each row, column, and box in the board
  """
  # iterate through ith row, column, and box
  for i in range(0, 9):
    # run hidden
    poss = hidden(poss, row_coordinates(i))
    poss = hidden(poss, col_coordinates(i))
    poss = hidden(poss, box_coordinates(i))
  return poss

@app.route("/", methods=["GET", "POST"])
def index():
  if request.method == "POST":
    global sudoku_board
    global poss
    global numsteps
    global coords
    global explanations
    global poss_h
    global something_changed
    global explanations_h
    global boxes

    # generate list of combinations for hidden singles, doubles, triples, quadruples
    combos = []
    # iterate through combination size
    for i in range(1, 5):
    # add all combinations
        combos += list(itertools.combinations(range(1, 10), i))
    
    boxes = [[n, m] for n in range(0, 9, 3) for m in range(0, 9, 3)]
    poss = [[[1, 2, 3, 4, 5, 6, 7, 8, 9] for j in range(9)] for i in range(9)]
    sudoku_board = []

    # Define sudoku_board[n][m] as the user-input value of the input with name nm 
    # (a convention we established in the html), if the value input is a number. 
    # If it is blank or non-numeric, set it to 0.
    for row in range(0,9):
        sudoku_board.append([])
        for column in range(0,9):
            x = request.form.get(str(row)+str(column))
            if x.isnumeric():
                sudoku_board[row].append(int(x))
            else:
                sudoku_board[row].append(0)
    # For each non-zero value in sudoku_board, set the corresponding possibility set to that value only
    for row in range(0, 9):
        for column in range(0, 9):
            p = sudoku_board[row][column]
            if p != 0:
                poss[row][column] = [p]

    #poss_h is the history of all poss arrays. highlights_ h is the history of all highlights arrays. etc.
    poss_h = []
    explanations_h = []

    for turn in range(100):
        # run all the functions
        something_changed = False
        explanations = [[[] for j in range(9)] for i in range(9)]
        poss = naked_single(poss)
        poss = naked_pair(poss)
        coords = coordinates(poss)
        poss = run_hidden(poss)
        poss = run_pointing(poss)
        poss_copy = copy.deepcopy(poss)
        poss_h.append(poss_copy)
        explanations_h.append(explanations)

        # check if nothing is changing. If so, we have reached the final step
        if not something_changed:
          numsteps = turn
          break
        
        # check if the sudoku is unsolvable,
        # which occurs if and only if a cell has no remaining possibilities
        for row in range(9):
          for col in range(9):
            if len(poss[row][col]) == 0:
              return render_template("apology.html")
    return render_template("output.html", step = 0, numsteps = numsteps, poss = poss_h[0], explanations = explanations_h[0])
  else:
      return render_template("index.html")
      


@app.route("/output", methods=["GET"])
def output():
  if request.method == "GET":
    # Figure out what step the user wants via GET, then navigate to that step by passing it to output.html
    step = int(request.args.get("step"))
    poss = poss_h[step]
    explanations = explanations_h[step]
    return render_template("output.html", step = step, numsteps = numsteps, poss = poss, explanations = explanations)