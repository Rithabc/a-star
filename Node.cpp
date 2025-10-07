
#include <stdio.h>
#include <iostream>
#include <vector>
#include <math.h>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <random>
#include <unordered_map>
#include <queue>
#include <stack>
#include <set>
#define INT_MAX = 2 ^ 9;
#define FLT_MAX 3.402823422e+38F
#include <cstring>

// Hash function for pair<int, int>
namespace std {
    template <>
    struct hash<std::pair<int, int>> {
        std::size_t operator()(const std::pair<int, int>& p) const noexcept {
            return std::hash<int>()(p.first) ^ (std::hash<int>()(p.second) << 1);
        }
    };
}

using namespace std;

std::mt19937 gen(std::chrono::system_clock::now().time_since_epoch().count());

int generateRandom(int low, int high)
{

    std::uniform_int_distribution<> dist(low, high);
    return dist(gen);
}

void dfs(int row, int col, int n, vector<vector<int>> &grid)
{
    int z = 0;
    if (row < 0 || row >= n || col < 0 || col >= n)
        return;
    if (grid[row][col] == 0)
        return;
    if (row - 1 >= 0 && row - 1 < n)
    {
        if (grid[row - 1][col] == 0)
        {
            z++;
            if (z > 1)
                return;
        }
    }
    if (row + 1 >= 0 && row + 1 < n)
    {
        if (grid[row + 1][col] == 0)
        {
            z++;
            if (z > 1)
                return;
        }
    }
    if (col + 1 >= 0 && col + 1 < n)
    {
        if (grid[row][col + 1] == 0)
        {
            z++;
            if (z > 1)
                return;
        }
    }
    if (col - 1 >= 0 && col - 1 < n)
    {
        if (grid[row][col - 1] == 0)
        {
            z++;
            if (z > 1)
                return;
        }
    }

    if (z > 1)
        return;
    grid[row][col] = 0;

    int cell = generateRandom(0, 3);
    switch (cell)
    {
    case 0:
    {
        dfs(row + 1, col, n, grid);
        dfs(row - 1, col, n, grid);
        dfs(row, col - 1, n, grid);
        dfs(row, col + 1, n, grid);
        break;
    }
    case 1:
    {
        dfs(row - 1, col, n, grid);
        dfs(row + 1, col, n, grid);
        dfs(row, col + 1, n, grid);
        dfs(row, col - 1, n, grid);
        break;
    }
    case 2:
    {
        dfs(row, col - 1, n, grid);
        dfs(row, col + 1, n, grid);
        dfs(row + 1, col, n, grid);
        dfs(row - 1, col, n, grid);
        break;
    }
    case 3:
    {
        dfs(row, col + 1, n, grid);
        dfs(row, col - 1, n, grid);
        dfs(row + 1, col, n, grid);
        dfs(row - 1, col, n, grid);
        break;
    }
    }
}

vector<vector<int>> generateGrid(int n)
{
    vector<vector<int>> grid(n, vector<int>(n, 1));
    int ind = generateRandom(0, n - 1);
    cout << ind << endl;
    // grid[ind][ind] = 0;
    dfs(ind, ind, n, grid);

    for (int row = 0; row < n; row++)
    {
        for (int col = 0; col < n; col++)
        {
            cout << grid[row][col] << " ";
        }
        cout << endl;
    }
    cout << "Map generated" <<endl;
    return grid;
}
// A*star professor's stratergy

int manhattan(int r1, int c1, int r2, int c2)
{
    return abs(r1 - r2) + abs(c1 - c2);
}

// A C++ Program to implement A* Search Algorithm

#define ROW 50
#define COL 50

// int ROW = 2;
// int COL = 2;

// Creating a shortcut for int, int pair type
typedef pair<int, int> Pair;

// Creating a shortcut for pair<int, pair<int, int>> type
typedef pair<double, pair<int, int>> pPair;

// A structure to hold the necessary parameters
struct cell
{
    // Row and Column index of its parent
    // Note that 0 <= i <= ROW-1 & 0 <= j <= COL-1
    int parent_i, parent_j;
    // f = g + h
    double f, g, h;
};

// A Utility Function to check whether given cell (row, col)
// is a valid cell or not.
bool isValid(int row, int col)
{
    // Returns true if row number and column number
    // is in range
    return (row >= 0) && (row < ROW) && (col >= 0) && (col < COL);
}

// A Utility Function to check whether the given cell is
// blocked or not
bool isUnBlocked(int grid[][COL], int row, int col)
{
    // Returns true if the cell is not blocked else false
    if (grid[row][col] == 0)
        return (true);
    else
        return (false);
}

// A Utility Function to check whether destination cell has
// been reached or not
bool isDestination(int row, int col, Pair dest)
{
    if (row == dest.first && col == dest.second)
        return (true);
    else
        return (false);
}

// A Utility Function to calculate the 'h' heuristics.
double calculateHValue(int row, int col, Pair dest)
{
    // Return using the distance formula
    return ((double)sqrt(
        (row - dest.first) * (row - dest.first) + (col - dest.second) * (col - dest.second)));
}

// A Utility Function to trace the path from the source
// to destination

vector<Pair> edit;
vector<string> directions;
void tracePath(cell cellDetails[][COL], Pair dest)
{
    // printf("\nThe Path is ");
    int row = dest.first;
    int col = dest.second;

    stack<Pair> Path;

    while (!(cellDetails[row][col].parent_i == row && cellDetails[row][col].parent_j == col))
    {
        Path.push(make_pair(row, col));
        int temp_row = cellDetails[row][col].parent_i;
        int temp_col = cellDetails[row][col].parent_j;
        row = temp_row;
        col = temp_col;
    }

    int tr = row;
    int tc = col;

    Path.push(make_pair(row, col));
    while (!Path.empty())
    {
        pair<int, int> p = Path.top();
        Path.pop();
        // cout << p.first << "," << p.second << " ";
        // cout << " " << tr << " " << tc << endl;
        if (p.first == tr - 1 && p.second == tc)
        {
            // printf("-> UP ");
            edit.push_back(make_pair(-1, 0));
            directions.push_back("UP");
        }
        else if (p.first == tr + 1 && p.second == tc)
        {
            // printf("-> DOWN ");
            edit.push_back(make_pair(1, 0));
            directions.push_back("DOWN");
        }
        else if (p.first == tr && p.second == tc - 1)
        {
            // printf("-> LEFT ");
            edit.push_back(make_pair(0, -1));
            directions.push_back("LEFT");
        }
        else if (p.first == tr && p.second == tc + 1)
        {
            // printf("-> RIGHT ");
            edit.push_back(make_pair(0, 1));
            directions.push_back("RIGHT");
        }
        tr = p.first;
        tc = p.second;
        // cout<<endl;

        // printf("-> (%d,%d) ", p.first, p.second);
    }

    return;
}

// A Function to find the shortest path between
// a given source cell to a destination cell according
// to A* Search Algorithm
void aStarSearch(int grid[][COL], Pair src, Pair dest)
{
    // If the source is out of range
    if (isValid(src.first, src.second) == false)
    {
        printf("Source is invalid\n");
        return;
    }

    // If the destination is out of range
    if (isValid(dest.first, dest.second) == false)
    {
        printf("Destination is invalid\n");
        return;
    }

    // Either the source or the destination is blocked
    if (isUnBlocked(grid, src.first, src.second) == false || isUnBlocked(grid, dest.first, dest.second) == false)
    {
        printf("Source or the destination is blocked\n");
        return;
    }

    // If the destination cell is the same as source cell
    if (isDestination(src.first, src.second, dest) == true)
    {
        // printf("We are already at the destination\n");
        return;
    }

    // Create a closed list and initialise it to false which
    // means that no cell has been included yet This closed
    // list is implemented as a boolean 2D array
    bool closedList[ROW][COL];
    memset(closedList, false, sizeof(closedList));

    // Declare a 2D array of structure to hold the details
    // of that cell
    cell cellDetails[ROW][COL];

    int i, j;

    for (i = 0; i < ROW; i++)
    {
        for (j = 0; j < COL; j++)
        {
            cellDetails[i][j].f = FLT_MAX;
            cellDetails[i][j].g = FLT_MAX;
            cellDetails[i][j].h = FLT_MAX;
            cellDetails[i][j].parent_i = -1;
            cellDetails[i][j].parent_j = -1;
        }
    }

    // Initialising the parameters of the starting node
    i = src.first, j = src.second;
    cellDetails[i][j].f = 0.0;
    cellDetails[i][j].g = 0.0;
    cellDetails[i][j].h = 0.0;
    cellDetails[i][j].parent_i = i;
    cellDetails[i][j].parent_j = j;

    /*
     Create an open list having information as-
     <f, <i, j>>
     where f = g + h,
     and i, j are the row and column index of that cell
     Note that 0 <= i <= ROW-1 & 0 <= j <= COL-1
     This open list is implemented as a set of pair of
     pair.*/
    set<pPair> openList;

    // Put the starting cell on the open list and set its
    // 'f' as 0
    openList.insert(make_pair(0.0, make_pair(i, j)));

    // We set this boolean value as false as initially
    // the destination is not reached.
    bool foundDest = false;

    while (!openList.empty())
    {
        pPair p = *openList.begin();

        // Remove this vertex from the open list
        openList.erase(openList.begin());

        // Add this vertex to the closed list
        i = p.second.first;
        j = p.second.second;
        closedList[i][j] = true;

        /*
         Generating all the 8 successor of this cell

             N.W   N   N.E
               \   |   /
                \  |  /
             W----Cell----E
                  / | \
                /   |  \
             S.W    S   S.E

         Cell-->Popped Cell (i, j)
         N -->  North       (i-1, j)
         S -->  South       (i+1, j)
         E -->  East        (i, j+1)
         W -->  West           (i, j-1)
         N.E--> North-East  (i-1, j+1)
         N.W--> North-West  (i-1, j-1)
         S.E--> South-East  (i+1, j+1)
         S.W--> South-West  (i+1, j-1)*/

        // To store the 'g', 'h' and 'f' of the 8 successors
        double gNew, hNew, fNew;

        //----------- 1st Successor (North) ------------

        // Only process this cell if this is a valid one
        if (isValid(i - 1, j) == true)
        {
            // If the destination cell is the same as the
            // current successor
            if (isDestination(i - 1, j, dest) == true)
            {
                // Set the Parent of the destination cell
                cellDetails[i - 1][j].parent_i = i;
                cellDetails[i - 1][j].parent_j = j;
                // printf("The destination cell is found\n");
                tracePath(cellDetails, dest);
                foundDest = true;
                return;
            }
            // If the successor is already on the closed
            // list or if it is blocked, then ignore it.
            // Else do the following
            else if (closedList[i - 1][j] == false && isUnBlocked(grid, i - 1, j) == true)
            {
                gNew = cellDetails[i][j].g + 1.0;
                hNew = calculateHValue(i - 1, j, dest);
                fNew = gNew + hNew;

                // If it isn’t on the open list, add it to
                // the open list. Make the current square
                // the parent of this square. Record the
                // f, g, and h costs of the square cell
                //                OR
                // If it is on the open list already, check
                // to see if this path to that square is
                // better, using 'f' cost as the measure.
                if (cellDetails[i - 1][j].f == FLT_MAX || cellDetails[i - 1][j].f > fNew)
                {
                    openList.insert(make_pair(
                        fNew, make_pair(i - 1, j)));

                    // Update the details of this cell
                    cellDetails[i - 1][j].f = fNew;
                    cellDetails[i - 1][j].g = gNew;
                    cellDetails[i - 1][j].h = hNew;
                    cellDetails[i - 1][j].parent_i = i;
                    cellDetails[i - 1][j].parent_j = j;
                }
            }
        }

        //----------- 2nd Successor (South) ------------

        // Only process this cell if this is a valid one
        if (isValid(i + 1, j) == true)
        {
            // If the destination cell is the same as the
            // current successor
            if (isDestination(i + 1, j, dest) == true)
            {
                // Set the Parent of the destination cell
                cellDetails[i + 1][j].parent_i = i;
                cellDetails[i + 1][j].parent_j = j;
                // printf("The destination cell is found\n");
                tracePath(cellDetails, dest);
                foundDest = true;
                return;
            }
            // If the successor is already on the closed
            // list or if it is blocked, then ignore it.
            // Else do the following
            else if (closedList[i + 1][j] == false && isUnBlocked(grid, i + 1, j) == true)
            {
                gNew = cellDetails[i][j].g + 1.0;
                hNew = calculateHValue(i + 1, j, dest);
                fNew = gNew + hNew;

                // If it isn’t on the open list, add it to
                // the open list. Make the current square
                // the parent of this square. Record the
                // f, g, and h costs of the square cell
                //                OR
                // If it is on the open list already, check
                // to see if this path to that square is
                // better, using 'f' cost as the measure.
                if (cellDetails[i + 1][j].f == FLT_MAX || cellDetails[i + 1][j].f > fNew)
                {
                    openList.insert(make_pair(
                        fNew, make_pair(i + 1, j)));
                    // Update the details of this cell
                    cellDetails[i + 1][j].f = fNew;
                    cellDetails[i + 1][j].g = gNew;
                    cellDetails[i + 1][j].h = hNew;
                    cellDetails[i + 1][j].parent_i = i;
                    cellDetails[i + 1][j].parent_j = j;
                }
            }
        }

        //----------- 3rd Successor (East) ------------

        // Only process this cell if this is a valid one
        if (isValid(i, j + 1) == true)
        {
            // If the destination cell is the same as the
            // current successor
            if (isDestination(i, j + 1, dest) == true)
            {
                // Set the Parent of the destination cell
                cellDetails[i][j + 1].parent_i = i;
                cellDetails[i][j + 1].parent_j = j;
                // printf("The destination cell is found\n");
                tracePath(cellDetails, dest);
                foundDest = true;
                return;
            }

            // If the successor is already on the closed
            // list or if it is blocked, then ignore it.
            // Else do the following
            else if (closedList[i][j + 1] == false && isUnBlocked(grid, i, j + 1) == true)
            {
                gNew = cellDetails[i][j].g + 1.0;
                hNew = calculateHValue(i, j + 1, dest);
                fNew = gNew + hNew;

                // If it isn’t on the open list, add it to
                // the open list. Make the current square
                // the parent of this square. Record the
                // f, g, and h costs of the square cell
                //                OR
                // If it is on the open list already, check
                // to see if this path to that square is
                // better, using 'f' cost as the measure.
                if (cellDetails[i][j + 1].f == FLT_MAX || cellDetails[i][j + 1].f > fNew)
                {
                    openList.insert(make_pair(
                        fNew, make_pair(i, j + 1)));

                    // Update the details of this cell
                    cellDetails[i][j + 1].f = fNew;
                    cellDetails[i][j + 1].g = gNew;
                    cellDetails[i][j + 1].h = hNew;
                    cellDetails[i][j + 1].parent_i = i;
                    cellDetails[i][j + 1].parent_j = j;
                }
            }
        }

        //----------- 4th Successor (West) ------------

        // Only process this cell if this is a valid one
        if (isValid(i, j - 1) == true)
        {
            // If the destination cell is the same as the
            // current successor
            if (isDestination(i, j - 1, dest) == true)
            {
                // Set the Parent of the destination cell
                cellDetails[i][j - 1].parent_i = i;
                cellDetails[i][j - 1].parent_j = j;
                // printf("The destination cell is found\n");
                tracePath(cellDetails, dest);
                foundDest = true;
                return;
            }

            // If the successor is already on the closed
            // list or if it is blocked, then ignore it.
            // Else do the following
            else if (closedList[i][j - 1] == false && isUnBlocked(grid, i, j - 1) == true)
            {
                gNew = cellDetails[i][j].g + 1.0;
                hNew = calculateHValue(i, j - 1, dest);
                fNew = gNew + hNew;

                // If it isn’t on the open list, add it to
                // the open list. Make the current square
                // the parent of this square. Record the
                // f, g, and h costs of the square cell
                //                OR
                // If it is on the open list already, check
                // to see if this path to that square is
                // better, using 'f' cost as the measure.
                if (cellDetails[i][j - 1].f == FLT_MAX || cellDetails[i][j - 1].f > fNew)
                {
                    openList.insert(make_pair(
                        fNew, make_pair(i, j - 1)));

                    // Update the details of this cell
                    cellDetails[i][j - 1].f = fNew;
                    cellDetails[i][j - 1].g = gNew;
                    cellDetails[i][j - 1].h = hNew;
                    cellDetails[i][j - 1].parent_i = i;
                    cellDetails[i][j - 1].parent_j = j;
                }
            }
        }
    }

    // When the destination cell is not found and the open
    // list is empty, then we conclude that we failed to
    // reach the destination cell. This may happen when the
    // there is no way to destination cell (due to
    // blockages)
    if (foundDest == false)
        printf("Failed to find the Destination Cell\n");

    return;
}

int minVal = 9999999;

int algoUtil(int grid[][COL], Pair src, Pair dest, vector<pair<int, int>> bots, stack<Pair> temp)
{
    directions.clear();
    while (!temp.empty())
    {
        // cout<<temp[0].first<<" "<<temp[0].second<<endl;
        // int random = generateRandom(0, temp.size() - 1);
        aStarSearch(grid, temp.top(), dest);
        if(directions.size()>minVal){
            directions.clear();
            edit.clear();
            return 9999999;
        }
        
        for (int i = 0; i < bots.size(); i++)
        {
            for (int j = 0; j < edit.size(); j++)
            {
                if (bots[i].first == dest.first && bots[i].second == dest.second)
                {
                }
                else
                {
                    if (bots[i].first + edit[j].first < 0 || bots[i].first + edit[j].first >= ROW || bots[i].second + edit[j].second < 0 || bots[i].second + edit[j].second >= COL)
                        continue;

                    bots[i].first += edit[j].first;
                    bots[i].second += edit[j].second;
                    if (grid[bots[i].first][bots[i].second] == 1)
                    {
                        bots[i].first -= edit[j].first;
                        bots[i].second -= edit[j].second;
                    }
                }
            }
        }
        // cout<<endl;
        temp = stack<Pair>();
        for (int i = 0; i < bots.size(); i++)
        {
            if (bots[i].first == dest.first && bots[i].second == dest.second)
            {
            }
            else
            {
                // cout<<"("<<bots[i].first<<","<<bots[i].second<<") ";
                temp.push(bots[i]);
            }
        }
        edit.clear();

        // for(int i = 0; i<temp.size(); i++){
        //     cout<<"("<<temp[i].first<<","<<temp[i].second<<") ";
        // }
    }
    // for(int i=0;i<directions.size();i++){
    //     cout<<directions[i]<<" ";
    // }
    // cout<<endl;
    // cout<<directions.size()<<endl;
    cout << src.first << "," << src.second << " to " << dest.first << "," << dest.second << " : ";
    cout << directions.size() << endl;
    if (directions.size() < minVal)
        minVal = directions.size();
    return directions.size();
}

struct hash_pair
{
    /* data */
    Pair p;
    vector<string> v;
    int len;
    hash_pair(Pair p, vector<string> v, int len) : p(p), v(v), len(len){};
};

void mapDfs(int i,int j,int prevI,int prevJ,int grid[][COL],vector<string> &ansDfs,vector<vector<int>> &visDfs,vector<Pair> &bots,unordered_map<pair<int,int>,int> &mpp,bool track = false){
    if(bots.size() == 0) return;
    if(mpp.size() == 0) return;
    if(i<0||j<0||i>=ROW||j>=COL) return;
    if(visDfs[i][j] == 1) return;
    visDfs[i][j]++;
    if(grid[i][j] == 1) return;

    // cout << i << "," << j << "->";

    if(!track){
    mpp.erase(make_pair(i,j));
    
    }
    // if( !bots.empty() && (i-1<0 || grid[i-1][j] == 1 || visDfs[i-1][j] == 1) && (i+1>ROW || grid[i+1][j] == 1 || visDfs[i+1][j] == 1) && (j-1<0 || grid[i][j-1] == 1 || visDfs[i][j-1] == 1) && (j+1>COL || grid[i][j+1] == 1 || visDfs[i][j+1] == 1) ){
    //     cout<<endl;
    //     cout << "if"<< endl;
    //     visDfs[prevI][prevJ] = 0;
    //     mapDfs(prevI,prevJ,i,j,grid,ansDfs,visDfs,bots);
    // }

    if(track){
        // cout << "backtrack " <<mpp.size() <<endl;
        if(prevI == i-1) ansDfs.push_back("DOWN");
        if(prevI == i+1) ansDfs.push_back("UP");
        if(prevJ == j-1) ansDfs.push_back("RIGHT");
        if(prevJ == j+1) ansDfs.push_back("LEFT");
        return;
    }

    if(prevI == i-1) ansDfs.push_back("DOWN");
    if(prevI == i+1) ansDfs.push_back("UP");
    if(prevJ == j-1) ansDfs.push_back("RIGHT");
    if(prevJ == j+1) ansDfs.push_back("LEFT");


    
    

    if(j+1<COL && visDfs[i][j+1] == 0 && grid[i][j+1] ==0){
        // ansDfs.push_back("RIGHT");
        mapDfs(i,j+1,i,j,grid,ansDfs,visDfs,bots,mpp);
    }
    if(i-1>=0 && visDfs[i-1][j] == 0 && grid[i-1][j]==0){
        // ansDfs.push_back("UP");
        mapDfs(i-1,j,i,j,grid,ansDfs,visDfs,bots,mpp);
    }
    if(j-1>=0 && visDfs[i][j-1] == 0 && grid[i][j-1] ==0){
        // ansDfs.push_back("LEFT");
        mapDfs(i,j-1,i,j,grid,ansDfs,visDfs,bots,mpp);
    }
    if(i+1<ROW && visDfs[i+1][j] == 0 && grid[i+1][j]==0){
        // ansDfs.push_back("DOWN");
        mapDfs(i+1,j,i,j,grid,ansDfs,visDfs,bots,mpp);
    }
    // mapDfs(i,j+1,i,j,grid,ansDfs,visDfs,bots);
    // mapDfs(i-1,j,i,j,grid,ansDfs,visDfs,bots);
    // mapDfs(i,j-1,i,j,grid,ansDfs,visDfs,bots);
    // mapDfs(i+1,j,i,j,grid,ansDfs,visDfs,bots);
    
    if(mpp.size() == 0) return;
    visDfs[prevI][prevJ] =0;
    mapDfs(prevI,prevJ,i,j,grid,ansDfs,visDfs,bots,mpp,true);
    
}

// Driver program to test above function
int main()
{
    /* Description of the Grid-
     1--> The cell is not blocked
     0--> The cell is blocked    */

    int grid[ROW][COL] ;

    vector<vector<int>> map = generateGrid(ROW);

    for(int i=0;i<ROW;i++){
        for(int j=0;j<COL;j++){
            grid[i][j] = map[i][j];
        }
    }
    
    vector<Pair> bots;
    
    unordered_map<pair<int,int>,int> mpp;



    for(int i=0;i<ROW;i++){
        for(int j=0;j<COL;j++){
            if(grid[i][j]==0){
            bots.push_back(make_pair(i,j));
            mpp[make_pair(i,j)]++;
            }
        }
    }
    vector<vector<int>> vis(ROW,vector<int>(ROW,0));
    vector<string> ans;
    queue<Pair> p;
    p.push(bots[0]);
    Pair dest;

    
    // // bots.clear();
    // while(!p.empty()){
    //     Pair t = p.front();
    //     p.pop();
    //     vis[t.first][t.second] = 1;
    //     //up 
    //     dest = t;
    //     int i = t.first;
    //     int j = t.second;

    //     if(j-1>=0 && vis[i][j-1] == 0 && grid[i][j-1] ==0){
    //         ans.push_back("RIGHT");
    //         p.push(make_pair(i,j-1));
    //     }

    //     if(j+1<COL && vis[i][j+1] == 0 && grid[i][j+1] == 0){
    //         ans.push_back("LEFT");
    //         p.push(make_pair(i,j+1));
    //     }

        
    //     if( i -1 >= 0 && vis[i-1][j] == 0 && grid[i-1][j]==0){
    //         ans.push_back("UP");
    //         p.push(make_pair(i-1,j));
    //     }

    //     if( i+1<ROW && vis[i+1][j] == 0 && grid[i+1][j]==0){
    //         ans.push_back("DOWN");
    //         p.push(make_pair(i+1,j));
    //     }

        
        
    // }
    // cout<<ans.size();
    cout<<endl;
    for(string s:ans){
        cout << s << " ";
    }
    cout<<endl;



    cout<< dest.first << "," << dest.second;

    cout << endl;
    
    vector<string> ansDfs;
    vector<vector<int>> visDfs(ROW,vector<int>(ROW,0));
    // mapDfs(bots[0].first,bots[0].second,-10,-10,grid,ansDfs,visDfs,bots,mpp);
    

    stack<Pair> st;
    st.push(bots[0]);
    int prevI = bots[0].first;
    int prevJ = bots[0].second;
    while(!mpp.size() == 0){
        Pair t = st.top();
        // cout << t.first<< "," <<t.second<<"->";
        mpp.erase(t);
        int i = t.first;
        int j = t.second;
        visDfs[i][j] = 1;

        if(prevI == i-1) ansDfs.push_back("DOWN");
        if(prevI == i+1) ansDfs.push_back("UP");
        if(prevJ == j-1) ansDfs.push_back("RIGHT");
        if(prevJ == j+1) ansDfs.push_back("LEFT");

        prevI = i;
        prevJ =j;

        if(i-1>=0 && visDfs[i-1][j] == 0 && grid[i-1][j]==0){
            st.push(make_pair(i-1,j));
        }
        else if(i+1<ROW && visDfs[i+1][j] == 0  && grid[i+1][j]==0){
            st.push(make_pair(i+1,j));
        }
        else if(j-1>=0 && visDfs[i][j-1] == 0  && grid[i][j-1]==0){
            st.push(make_pair(i,j-1));
        }
        else if(j+1<COL && visDfs[i][j+1] == 0 && grid[i][j+1]==0){
            st.push(make_pair(i,j+1));
        }else{
            // cout << i << ","<<j;
            // cout << "pop";
            st.pop();
        }

        

    }

    // for(string s:ansDfs){
    //     cout << s << " ";
    // }
    cout <<endl;
    cout << bots[1].first << "," << bots[1].second <<endl;
    cout << ansDfs.size();       

    
    
    return (0);
}