// Project Identifier: 5949F553E20B650AB0FB2266D3C0822B13D248B0
#include <iostream>
#include <string>
#include <algorithm>
#include <vector>
#include <limits>
#include <list>
#include <queue>
#include <cmath>
#include <iomanip>
#include <numeric>
#include "getopt.h"

using namespace std;

enum class Mode
{
    None = 0,
    SpanningTree,
    Fast,
    Optimal,
};

struct Pokemon
{
int xCoord;
int yCoord;
Pokemon(int x, int y);
};

Pokemon::Pokemon(int x, int y)
:xCoord(x),yCoord(y){}

struct MSTInfo
{
  bool visited = false;
  double minWeightSquared = numeric_limits<double>::infinity();
  int parentVertex = -1;
};

class Solver
{
  public:
  Solver(string mode);
  void readPokemon();
  void solve();

  private:
  vector<Pokemon> pokemon;
  vector<MSTInfo> MSTinformation;
  vector<int>fastTSPTour;
  vector<int> bestOPTTSPTour;
  std::vector<std::vector<double>> OPTTSPDistances;
  string mode;
  int numPokemon;
  double bestOPTTSPCost = numeric_limits<double>::infinity();
  bool verifyMSTConnection(const Pokemon& p1, const Pokemon& p2);
  double calculateDistance(const Pokemon& p1, const Pokemon& p2);
  double calculateSquaredDistance(const Pokemon& p1, const Pokemon& p2);
  void solveMST();
  double MSTCostForSubset(vector<int>& nodes);
  double insertionCost(int pokemonToInsert, int prevPokemon, int nextPokemon);
  double solveFASTTSP(bool printResult);
  void calcDistanceMatrix();
  double calcPartialPathCost(vector<int>&path, size_t permLength);
  double calcLowerBound(vector<int>& path, size_t permLength);
  bool promisingPath(vector<int>& path, size_t permLength);
  void processSolution(vector<int>& path);
  void genPerms(vector<int>& path, size_t permLength);
  void solveOPTTSP();
};

Solver::Solver(string mode)
:mode(mode){}

void Solver::readPokemon()
{
  cin >> numPokemon;
  pokemon.reserve(static_cast<size_t>(numPokemon));
  if(mode == "MST") MSTinformation.reserve(static_cast<size_t>(numPokemon));
  for(int i=0; i<numPokemon; i++)
  {
    int x,y;
    cin >> x >> y;
    Pokemon p(x,y);
    pokemon.push_back(p);
    if(mode == "MST") 
    {
      MSTInfo info;
      MSTinformation.push_back(info);
    }
  }
}

void Solver::solve()
{
  if(mode == "MST")
  {
    solveMST();
  }
  else if(mode == "FASTTSP")
  {
    solveFASTTSP(true);
  }
  else
  {
    solveOPTTSP();
  }
}

// double Solver::calculateSquaredDistance(const Pokemon& p1, const Pokemon& p2)
// {
//   long double xDiff = static_cast<long double>(p1.xCoord) - p2.xCoord;
//   long double yDiff = static_cast<long double>(p1.yCoord) - p2.yCoord;

//   return static_cast<double>(xDiff*xDiff + yDiff*yDiff);
// }

inline double Solver::calculateSquaredDistance(const Pokemon& p1, const Pokemon& p2)
{
  long double xDiff = static_cast<long double>(p1.xCoord) - p2.xCoord;
  long double yDiff = static_cast<long double>(p1.yCoord) - p2.yCoord;
  return static_cast<double>(xDiff*xDiff + yDiff*yDiff);
}

inline double Solver::calculateDistance(const Pokemon& p1, const Pokemon& p2)
{
  long double xDiff = static_cast<long double>(p1.xCoord) - p2.xCoord;
  long double yDiff = static_cast<long double>(p1.yCoord) - p2.yCoord;

  return static_cast<double>(sqrt(xDiff*xDiff + yDiff*yDiff));
}

bool Solver::verifyMSTConnection(const Pokemon& p1, const Pokemon& p2)
{
    // Check if both are in the "sea" (both x and y are negative for each Pokemon)
    bool p1IsSea = (p1.xCoord < 0 && p1.yCoord < 0);
    bool p2IsSea = (p2.xCoord < 0 && p2.yCoord < 0);
    
    if (p1IsSea && p2IsSea) {
        return true;
    }

    // Check if both are on "land" (either x or y is positive for each Pokemon)
    bool p1IsLand = (p1.xCoord > 0 || p1.yCoord > 0);
    bool p2IsLand = (p2.xCoord > 0 || p2.yCoord > 0);

    if (p1IsLand && p2IsLand) {
        return true;
    }

    // Check if one of them is on the "coast"
    bool p1IsCoast = (p1.yCoord == 0 && p1.xCoord <= 0) || (p1.xCoord == 0 && p1.yCoord <= 0);
    bool p2IsCoast = (p2.yCoord == 0 && p2.xCoord <= 0) || (p2.xCoord == 0 && p2.yCoord <= 0);

    if (p1IsCoast || p2IsCoast) {
        return true;
    }

    // If none of the above conditions match, the connection is invalid
    return false;
}

void Solver::solveMST()
{
  auto np = static_cast<size_t>(numPokemon);
  double totalWeight = 0.0;
  MSTinformation[0].minWeightSquared = 0.0;
  for(size_t i=0; i<np; i++)
  {
    size_t smallestUnvisitedIndex = np;
    double smallestWeightSquared = numeric_limits<double>::infinity();
    //find smallest unvisited index
    for(size_t j=0; j<np; j++)
    {
      if(!MSTinformation[j].visited && 
      MSTinformation[j].minWeightSquared < smallestWeightSquared)
      {
        smallestUnvisitedIndex = j;
        smallestWeightSquared = MSTinformation[j].minWeightSquared;
      }
    }
    MSTinformation[smallestUnvisitedIndex].visited = true;
    if(smallestUnvisitedIndex != 0)
    {
      size_t parentIndex = static_cast<size_t>(MSTinformation[smallestUnvisitedIndex].parentVertex);
      if(verifyMSTConnection(pokemon[smallestUnvisitedIndex],
      pokemon[parentIndex]))
      {
        //add this node to mst if the connection is valid
        double distanceSquared = calculateSquaredDistance(pokemon[smallestUnvisitedIndex],
        pokemon[parentIndex]);
        totalWeight += sqrt(distanceSquared);
      }
    }
    //recalculate the distance between new node and every unvisited node
    for(size_t j = 0; j < np; j++)
    {
      if(!MSTinformation[j].visited)
      {
        double distanceSquared = calculateSquaredDistance(pokemon[static_cast<size_t>(smallestUnvisitedIndex)], 
        pokemon[j]);
        if(verifyMSTConnection(pokemon[static_cast<size_t>(smallestUnvisitedIndex)],
        pokemon[j])
        && distanceSquared < MSTinformation[static_cast<size_t>(j)].minWeightSquared)
        {
          MSTinformation[j].minWeightSquared = distanceSquared;
          MSTinformation[j].parentVertex = static_cast<int>(smallestUnvisitedIndex);
        }
      }
    }
  }
  cout << totalWeight << "\n";
  for(int i = 1; i < numPokemon; i++)
  {
    int currentParentVertex = MSTinformation[(size_t)i].parentVertex;
    cout << min(i,currentParentVertex) << " " << max(i,currentParentVertex) << "\n";
  }
}


double Solver::insertionCost(int pokemonToInsert, int prevPokemon, int nextPokemon)
{
  size_t pti = static_cast<size_t>(pokemonToInsert);
  size_t pp = static_cast<size_t>(prevPokemon);
  size_t np = static_cast<size_t>(nextPokemon);
  return (calculateDistance(pokemon[pp], pokemon[pti]) + 
  calculateDistance(pokemon[pti], pokemon[np]) - 
  calculateDistance(pokemon[pp], pokemon[np]));  
}

double Solver::solveFASTTSP(bool printResult)
{
  std::list<int> pokemonTour = {0,1,2};
  for(int pokemonToInsert = 3; pokemonToInsert < numPokemon; pokemonToInsert ++)
  {
    double minCost = std::numeric_limits<double>::infinity();
    std::list<int>::iterator bestInsertPosition;

    for(auto it = pokemonTour.begin(); it != pokemonTour.end(); ++it)
    {
      auto nextPos = std::next(it);
      if(nextPos == pokemonTour.end()) nextPos = pokemonTour.begin();
      double currentCost = insertionCost(pokemonToInsert,*it,*nextPos);
      if(currentCost < minCost)
      {
        minCost = currentCost;
        bestInsertPosition = nextPos;
      }
    }
    pokemonTour.insert(bestInsertPosition,pokemonToInsert);
  }
  double totalWeight = 0.0;
  auto currentPos = pokemonTour.begin();
  auto tourEnd = pokemonTour.end();
  while(true)
  {
    auto nextPos = std::next(currentPos);
    if(nextPos == tourEnd) nextPos = pokemonTour.begin();
    auto cp = static_cast<size_t>(*currentPos);
    auto np = static_cast<size_t>(*nextPos);
    totalWeight +=calculateDistance(pokemon[cp],pokemon[np]);
    ++currentPos;
    if(currentPos == tourEnd) break;
  }
  fastTSPTour = vector<int>(pokemonTour.begin(),pokemonTour.end());
  if(printResult)
  {
  cout << totalWeight << "\n";
  auto zeroIt = std::find(fastTSPTour.begin(), fastTSPTour.end(), 0);
  std::rotate(fastTSPTour.begin(), zeroIt, fastTSPTour.end());
  for(size_t i=0; i< fastTSPTour.size(); i++)
  {
    cout << fastTSPTour[i] << " ";
  }
  }
  return totalWeight;
}

double Solver::MSTCostForSubset(vector<int>& nodes)
{
    size_t n = nodes.size();
    if (n == 0) return 0.0;
    vector<double> minWeight(n, numeric_limits<double>::infinity());
    minWeight[0] = 0.0;

    double mstCost = 0.0;
    for (size_t remainingNodes = n; remainingNodes > 0; remainingNodes--) {

        double minCost = numeric_limits<double>::infinity();
        size_t u = 0;
        for (size_t j = 0; j < remainingNodes; j++) {
            if (minWeight[j] < minCost) {
                minCost = minWeight[j];
                u = j;
            }
        }
        mstCost += sqrt(minCost);
        std::swap(nodes[u], nodes[remainingNodes - 1]);
        std::swap(minWeight[u], minWeight[remainingNodes - 1]);

        for (size_t v = 0; v < remainingNodes - 1; v++) {
            double cost = OPTTSPDistances[static_cast<size_t>(nodes[remainingNodes - 1])][static_cast<size_t>(nodes[v])];
            minWeight[v] = min(minWeight[v], cost);
        }
    }
    return mstCost;
}

void Solver::calcDistanceMatrix()
{
  auto np = static_cast<size_t>(numPokemon);
  OPTTSPDistances.resize(np, std::vector<double>(np, numeric_limits<double>::infinity()));
  for (size_t i = 0; i < np; i++) 
  {
    for (size_t j = 0; j < np; j++) 
    {
      if(OPTTSPDistances[i][j] != numeric_limits<double>::infinity()) continue;

      OPTTSPDistances[i][j] = OPTTSPDistances[j][i] = (i != j) ? calculateSquaredDistance(pokemon[i], pokemon[j]) : 0.0;
    }
  }
}

double Solver::calcPartialPathCost(vector<int>&path, size_t permLength)
 {
  double cost = 0.0;
  for(size_t i = 0; i < permLength - 1; i++)
  {
    cost += sqrt(OPTTSPDistances[static_cast<size_t>(path[i])][static_cast<size_t>(path[i+1])]);
  }
  return cost;
 }
 
 double Solver::calcLowerBound(std::vector<int>& path, size_t permLength)
 {
   double bound = 0.0;
    for (size_t i = 1; i < permLength; i++) {
        bound += sqrt(OPTTSPDistances[static_cast<size_t>(path[i - 1])][static_cast<size_t>(path[i])]);
    }

    // Add return-to-start cost if the path is complete
    if (permLength == path.size()) {
        bound += sqrt(OPTTSPDistances[static_cast<size_t>(path.back())][static_cast<size_t>(path[0])]);
        return bound;
    }

    // Check early if bound already exceeds the best cost
    if (bound >= bestOPTTSPCost) return bound;

    // Find unvisited nodes
    vector<int> unvisited;
    for (size_t i = permLength; i < OPTTSPDistances.size(); i++)
    {
    unvisited.push_back(static_cast<int>(path[i]));
    }

    // Calculate MST cost for unvisited nodes
    double mstCost = MSTCostForSubset(unvisited);
    bound += mstCost;

    // Check again if bound already exceeds the best cost
    if (bound >= bestOPTTSPCost) return bound;

    // Add minimum connection costs for entry and exit
    double entryCost = numeric_limits<double>::infinity();
    double exitCost = numeric_limits<double>::infinity();
    for (int node : unvisited) {
        entryCost = min(entryCost, OPTTSPDistances[static_cast<size_t>(path[permLength - 1])][static_cast<size_t>(node)]);
        exitCost = min(exitCost, OPTTSPDistances[static_cast<size_t>(node)][static_cast<size_t>(path[0])]);
    }
    bound += sqrt(entryCost) + sqrt(exitCost);

    return bound;
}


bool Solver::promisingPath(vector<int>& path, size_t permLength)
{
if(permLength < 2) return true;
return calcLowerBound(path,permLength) <= bestOPTTSPCost;
}

void Solver::processSolution(std::vector<int>& path) 
{
  double totalCost = 0.0;
        // Calculate cost of complete tour including return to start
        for (size_t i = 0; i < path.size(); i++) {
            size_t current = static_cast<size_t>(path[i]);
            size_t next = static_cast<size_t>(path[(i + 1) % path.size()]);
            totalCost += sqrt(OPTTSPDistances[current][next]);
        }

        if (totalCost < bestOPTTSPCost) {
            bestOPTTSPCost = totalCost;
            bestOPTTSPTour = path;
        }
    }

void Solver::genPerms(vector<int>& path, size_t permLength)
{
  if(permLength == path.size())
  {
   processSolution(path);
   return; 
  }
  if(!promisingPath(path,permLength)) return;
  for(size_t i = permLength; i < path.size(); i++)
  {
    swap(path[permLength], path[i]);
    genPerms(path,permLength+1);
    swap(path[permLength], path[i]);
  }
}

void Solver::solveOPTTSP()
{
  calcDistanceMatrix();
  bestOPTTSPCost = solveFASTTSP(false);
  auto np = static_cast<size_t>(numPokemon);
  vector<int> path(np);
  std::iota(path.begin(), path.end(), 0);
  bestOPTTSPTour = fastTSPTour;
  genPerms(path,1);
  cout << bestOPTTSPCost << "\n";
  if(bestOPTTSPTour[0] != 0)
  {
  auto zeroIt = std::find(bestOPTTSPTour.begin(), bestOPTTSPTour.end(), 0);
  std::rotate(bestOPTTSPTour.begin(), zeroIt, bestOPTTSPTour.end());
  }
    for (size_t i = 0; i < bestOPTTSPTour.size(); i++) {
        cout << bestOPTTSPTour[i] << " ";
    }
}


void printHelp(char *argv[]) {
  cout << "Usage: " << argv[0] << " [-m resize|reserve|nosize] | -h\n";
  cout << "This program is to help you learn command-line processing,\n";
  cout << "reading data into a vector, the difference between resize and\n";
  cout << "reserve and how to properly read until end-of-file." << endl;
}  // printHelp()

struct Options
{
Mode mode = Mode::None;
};

void setModes(int argc, char * argv[], Options &options) {
  // These are used with getopt_long()
  opterr = false; // Let us handle all error output for command line options
  int choice;
  int index = 0;
  option long_options[] = {
    {"help", no_argument,         nullptr, 'h' },
    {"mode", required_argument,   nullptr, 'm' },
    {nullptr,  0,                 nullptr, '\0'},
  };  // long_options[]


  while ((choice = getopt_long(argc, argv, "hm:", long_options, &index)) != -1) {
    switch (choice) {
      case 'h':
      printHelp(argv);
      exit(0);

      case 'm':
      {
        string mode{optarg};
        if(mode != "MST" && mode != "FASTTSP" && mode != "OPTTSP")
        {
            cerr << "Error: Invalid mode" <<endl;
            exit(1);
        }
        if(mode == "MST")
        {
            options.mode = Mode::SpanningTree;
        }
        else if(mode == "FASTTSP")
        {
            options.mode = Mode::Fast;
        }
        else if(mode == "OPTTSP")
        {
            options.mode = Mode::Optimal;
        }
        break;
      }
      
      default:
      cerr << "Error: Invalid command line option"<<endl;
      exit(1);


    }  // switch ..choice
    if(options.mode == Mode::None)
    {
        cerr << "Error: No mode specified" << endl;
        exit(1);
    }
  }  // while
}  // getMode()

int main(int argc, char* argv[])
{
ios_base::sync_with_stdio(false);

cout << std::setprecision(2); //Always show 2 decimal places
cout << std::fixed; //Disable scientific notation for large numbers
Options options;
setModes(argc,argv,options);
string mode;
if(options.mode == Mode::SpanningTree) mode = "MST";
else if(options.mode == Mode::Fast) mode = "FASTTSP";
else mode = "OPTTSP";
Solver solver(mode);
solver.readPokemon();
solver.solve();
return 0;
}

