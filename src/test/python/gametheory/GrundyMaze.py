
'''
following 25.3 Sprague–Grundy theorem
from Competitive Programmer’s Handbook
by Antti Laaksonen
'''
class GrundyMaze :

    '''
    find the moves available to pos in grid.
    grid is a list of lists of grid values.
    pos is a tuple of (row, col).

    the moves available are steps left or up until a barrier is reached or boundary of grid.
    '''
    def findMoves(self, grid: list , visited: set, pos : tuple) -> list :
       # limit offsets.  restricted to moving
       moves = []
       m = len(grid)
       n = len(grid[0])
       # moving upward from pos until hit barrier or end
       for r2 in range(pos[0]-1, -1, -1):
           c2 = pos[1]
           if (r2 < 0) : break
           # here for rules about not able to move to this square
           # e.g. if grid == -10, cannot occupy square
           if grid[r2][c2] == -10 : break
           if (r2,c2) in visited : continue
           moves.append((r2, c2))
       # moving left until reach barrier or boundary
       for c2 in range(pos[1] - 1, -1, -1):
           r2 = pos[0]
           if (c2 < 0): break
           # here for rules about not able to move to this square
           # e.g. if grid == -10, cannot occupy square
           if grid[r2][c2] == -10: break
           if (r2, c2) in visited: continue
           moves.append((r2, c2))
       return moves

    '''
    pos is a tuple of (row, col)
    grid is a list of lists of grid values.
    '''
    def grundyNumber(self, grid: list, visited: set, pos : tuple = (-100000, -100000)) -> int :
        if pos[0] == -100000 :
            pos = (len(grid)-1, len(grid[0])-1)
        # this can be edited for context.  assuming grid has -1 for fillable,
        # and 0 = end state, grundyNumber = 0
        # and -10 for unoccupiable
        #if grid[pos[0]][pos[1]] != -1:
        #    return grid[pos[0]][pos[1]];
        visited.add(pos)
        moves = self.findMoves(grid, visited, pos)
        s = set()
        for x in moves :
           if x in visited:
               y = grid[x[0]][x[1]]
           else :
               y = self.grundyNumber(grid, visited, x)
           s.add(y)
        # calc mex, that is smallest number not in set s
        ret = 0
        while ret in s :
            ret += 1
        # update grid:
        grid[pos[0]][pos[1]] = ret
        return ret

    '''
    find the moves available to pos in grid.
    grid is a list of lists of grid values.
    pos is a tuple of (row, col)
    '''
    def findMoves0(self, grid: list , visited: set, pos : tuple) -> list :
       # limit offsets.  restricted to moving
       #offsets = [[0,-1],[-1,0], [0,1], [1,0]]
       offsets = [[1,0], [-1,0], [0,1], [0,-1]]
       moves = []
       m = len(grid)
       n = len(grid[0])
       for offset in offsets:
           r2 = pos[0] + offset[0]
           c2 = pos[1] + offset[1]
           if (r2 < 0 or c2 < 0 or r2 == m or c2 == n) : continue
           # here for rules about not able to move to this square
           # e.g. if grid == -10, cannot occupy square
           if (grid[r2][c2] == -10 or (r2, c2) in visited) : continue
           moves.append((r2, c2))
       return moves

if __name__ == '__main__':
    m = 5
    n = 5
    grid = [[0 for i in range(n)] for i in range(m)]
    grid[0][2] = -10
    grid[1][0] = -10
    grid[1][4] = -10
    grid[2][2] = -10
    grid[3][0] = -10
    visited = set()
    gr = GrundyMaze()
    # starts from bottom right of grid by default
    gr.grundyNumber(grid, visited)
    print(f'grid={grid}')