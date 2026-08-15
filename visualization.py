"""
1. Make additional spot for - or | elements 
to show the connection between residues.

2. Choose boundaries for future figure, because negative coordinates can exist.
"""

SYMBOLS = {
    "H": "H", 
    "P": "P"
           }

#To vizualize ASCII figure for conformation.
def ascii_vizualization(sequence, coords):
    if len(sequence) != len(coords):
        raise ValueError(f"Different lengths: {len(sequence) and len(coords)}")


    #min and max for a window
    xs = [x for x, _ in coords]
    ys = [y for _, y in coords]
    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)

    #for bonds 
    width = 2 * (max_x - min_x) + 1
    height = 2 * (max_y - min_y) + 1

    #table formation of same size as future protein representation
    grid = [[" "] * width for _ in range(height)]

    #Coordinates -> cell position (column, row).
    def pos_to_cell(x,y):
        column = 2 * (x - min_x)
        row = 2 * (y - min_y)
        return column, row

    #Draw residue in cell
    for (x,y), residue in zip(coords, sequence):
        column, row = pos_to_cell(x,y)
        grid[row][column] = SYMBOLS[residue]

    #Draw bonds
    for (x1, y1), (x2, y2) in zip(coords, coords[1:]):
        column1, row1 = pos_to_cell(x1, y1)
        column2, row2 = pos_to_cell(x2, y2)
        mid_column, mid_row = (column1 + column2) // 2, (row1 + row2) // 2
        grid[mid_row][mid_column] = "-" if row1 == row2 else "|"

    #Combine all rows in one figure
    return "\n".join("".join(row) for row in reversed(grid)) #to print from top and not from bottom line
