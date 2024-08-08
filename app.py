import os

from flask import Flask, flash, jsonify, redirect, render_template, request, session

# Configure application
app = Flask(__name__)

# Ensure templates are auto-reloaded
app.config["TEMPLATES_AUTO_RELOAD"] = True


#check if an array contains two of the same number -- currently not being used. May not be useful, who knows
def isoverlap(array):
    for n in range (0, 7):
        for m in range (n+1,8):
            if (array[n] == array[m]):
                return True
    return False

#check overlap in rows columns and boxes
def reduce_basic(poss):
    for n in range (0,9):
        for m in range (0,9):
            if int(len(poss[n][m])) == 1:
                x = poss[n][m][0]
                    #go through x's row, and remove x from the possibilities of all the other cells
                for m2 in range (0,9):
                    if m2 != m:
                        poss[n][m2].remove(x)
                #go through x's column, and remove x from the possibilities of all the other cells
                for n2 in range (0,9):
                    if n2 != n:
                        poss[n2][m].remove(x)
                #go through x's box, and remove x from the possibilities of all the other cells
                for n2 in range (n - n%3, n - n%3 + 3):
                    for m2 in range (m - m%3, m - m%3 + 3):
                        if n2 != n and m2 != m:
                            poss[n2][m2].remove(x)

    return poss

@app.route("/", methods=["GET", "POST"])
def index():
    if request.method == "POST":
        
        # TODO: Add the user's entry into a variable called sudoku board
        # At the same time, I am populating an array called "possibilities" with the possibilities of each cell.
        sudoku_board = []
        poss = []
        for row in range(0,9):
            sudoku_board.append([])
            poss.append([])
            for column in range(0,9):
                x = request.form.get(str(row)+str(column))
                if x.isnumeric():
                    sudoku_board[row].append(int(x))
                    poss[row].append([int(x)])
                else:
                    sudoku_board[row].append(0)
                    poss[row].append([1,2,3,4,5,6,7,8,9])
        
        #print sudoku board
        for n in range (0,9): 
            print (sudoku_board[n])

        #calculate possibilities for each cell (we're going to want to do this a lot of times)
        poss = reduce_basic(poss)

        #print possibilities
        for n in range (0,9): 
            print (poss[n])

        return redirect("/")

    else:
        
        # TODO: Display the entries in the database on index.html
        return render_template("index.html")
        #return render_template("index.html", names = db.execute(SELECT name FROM birthdays;))


