package COMP6651_Project;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

public class GraphAnalyzer {
    static class Vertex {
        int id;
        double distanceToDestination = 0.0;
        double highestCostToNode = 0.0;
        List<Vertex> currentPathToNode = new ArrayList<>();
        int degree = 0;
        double x;
        double y;
        List<Vertex> neighbors = new ArrayList<>();
        List<Double> distances = new ArrayList<>();
        boolean visited = false;
        double distance = Double.MAX_VALUE;
        Vertex predecessor = null;

        Vertex(int id) {
            this.id = id;
        }

        void addNeighbor(Vertex neighbor) {

            this.neighbors.add(neighbor);
            this.distances.add(1.0);
        }
        void addNeighbor(Vertex neighbor, double distance){
            this.neighbors.add(neighbor);
            this.distances.add(distance);
        }
    }

    private List<Vertex> vertices;
    private List<List<Vertex>> allComponents = new ArrayList<>();
    private double LmaxDFS = 0;
    private double LmaxDijkstra = 0;
    private double LmaxAStar = 0;
    private double LmaxBeamSearch = 0;
    private int maxDegreeInLSP_DFS = 0;
    private int maxDegreeInLSP_Dijkstra = 0;
    public GraphAnalyzer(String filename) throws IOException {
        vertices = new ArrayList<>();
        readGraphFromCSV(filename);
        findAllConnectedComponents();
    }

    private void readGraphFromCSV(String filename) throws IOException {
        List<String> lines = Files.readAllLines(Paths.get(filename));

        // Parse the first line for the number of vertices;
        int lineCount = 0;
        String[] firstLineTokens = lines.get(lineCount).split(" ");
        while(firstLineTokens[0].contains("%")){
            lineCount++;
             firstLineTokens = lines.get(lineCount).split(" ");
        }
        int numVertices = Integer.parseInt(firstLineTokens[1]); // Number of vertices USE index 1, middle value is vertices for MTX and EDGE files.
        for (int i = 0; i < numVertices; i++) {
            vertices.add(new Vertex(i)); // Initialize vertices
        }

        // Parse subsequent lines for edges
        for (int i = lineCount+1; i < lines.size(); i++) {
            String[] tokens = lines.get(i).split(" ");
            int to;
            int from;
            if(tokens.length == 2){
            from = Integer.parseInt(tokens[0]) - 1; // Adjust for 0-based indexing
             to = Integer.parseInt(tokens[1]) - 1; // Adjust for 0-based indexing
                vertices.get(from).addNeighbor(vertices.get(to));
                vertices.get(to).addNeighbor(vertices.get(from)); // Assuming undirected graph
            }
            if(tokens.length == 6){
                from = Integer.parseInt(tokens[0]) - 1; // Adjust for 0-based indexing
                to = Integer.parseInt(tokens[3]) - 1; // Adjust for 0-based indexing
                vertices.get(from).x = Double.parseDouble(tokens[1]);
                vertices.get(from).y = Double.parseDouble(tokens[2]);
                vertices.get(to).x = Double.parseDouble(tokens[4]);
                vertices.get(to).y = Double.parseDouble(tokens[5]);
                vertices.get(from).addNeighbor(vertices.get(to), Double.parseDouble(tokens[2]));
                vertices.get(to).addNeighbor(vertices.get(from), Double.parseDouble(tokens[2]));
            }
        }
    }

    private void findAllConnectedComponents() {
        boolean[] visited = new boolean[vertices.size()];
        for (Vertex v : vertices) {
            if (!visited[v.id]) {
                List<Vertex> component = new ArrayList<>();
                dfsFindComponent(v, visited, component);
                allComponents.add(component);
            }
        }
    }

    private void dfsFindComponent(Vertex v, boolean[] visited, List<Vertex> component) {
        visited[v.id] = true;
        component.add(v);
        for (Vertex neighbor : v.neighbors) {
            if (!visited[neighbor.id]) {
                dfsFindComponent(neighbor, visited, component);
            }
        }
    }

    private List<Vertex> getLargestComponent() {
        List<Vertex> lcc = allComponents.stream().max(Comparator.comparingInt(List::size)).orElse(new ArrayList<>());
        return lcc;
    }
    
    private int calculateDeltaLCC() {
 	   List<Vertex> lcc = getLargestComponent();
 	   int maxDegree = 0;
 	   for (Vertex v : lcc) {
 	       maxDegree = Math.max(maxDegree, v.neighbors.size()); 
 	   }
 	   return maxDegree;
 	}
 
	 private int calculateVLCC() {
	     return getLargestComponent().size();
	 }
	
	 private double calculateAverageDegreeLCC() {
	     List<Vertex> lcc = getLargestComponent();
	     return lcc.stream().mapToInt(v -> v.neighbors.size()).average().orElse(0);
	 }

    private void dfsCalculateDegree(Vertex v, int depth, Set<Vertex> visited) {
        // Add the current vertex to the visited set at the beginning of the method
        if (!visited.add(v)) {
            // If adding v returns false, it means v was already in the set and thus already visited on this path
            return; // Skip this vertex to avoid loops
        }

        // Update maxDegreeInLSP_DFS and LmaxDFS based on the current vertex
        maxDegreeInLSP_DFS = Math.max(maxDegreeInLSP_DFS, v.neighbors.size());
        LmaxDFS = Math.max(LmaxDFS, depth); 



        // Recursively visit neighbors
        for (Vertex neighbor : v.neighbors) {
            dfsCalculateDegree(neighbor, depth + 1, visited); // Pass the same visited set to the recursive call
        }
    }
    
    public void runDFSForDeltaLCC() {
        List<Vertex> lcc = getLargestComponent();
        maxDegreeInLSP_DFS = 0; 
        LmaxDFS = 0;
        for (Vertex vertex : vertices) {
            vertex.visited = false; // Reset visited for all vertices
        }
        for (Vertex vertex : lcc) {
            Set<Vertex> visited = new HashSet<>(); // Initialize the visited set for this starting vertex
            dfsCalculateDegree(vertex, 1, visited); // Start depth from 1
        }
    }

    private void dijkstraCalculateDegree(Vertex start, List<Vertex> lcc) {
        PriorityQueue<Vertex> pq = new PriorityQueue<>(Comparator.comparingDouble(v -> v.distance));
        start.distance = 0;
        pq.offer(start);

        while (!pq.isEmpty()) {
            Vertex u = pq.poll();
            
            for (Vertex v : u.neighbors) {
                if (!lcc.contains(v) || v.distance <= u.distance + 1) continue;
                
                double weight = 1.0; // Assuming unity weight for simplicity
                v.distance = u.distance + weight;
                pq.offer(v);
            }
        }
    }
    
    public void runAlgorithmsAndCalculateMetrics() {
      runDFSForDeltaLCC();
      runDijkstraForDeltaLCC();
  }
      
    public void runDijkstraForDeltaLCC() {
      List<Vertex> lcc = getLargestComponent();
      maxDegreeInLSP_Dijkstra = 0;
      LmaxDijkstra = 0; // Ensure this is reset only once before all runs

      // Loop through each vertex in the LCC as the starting point for Dijkstra's algorithm
      for (Vertex startVertex : lcc) {
          // Reset distances for all vertices before each run
          for (Vertex vertex : vertices) {
              vertex.distance = Double.MAX_VALUE;
              vertex.predecessor = null;
          }
          
          // Run Dijkstra from the current start vertex
          dijkstraCalculateDegree(startVertex, lcc);

          // After the run, check all vertices to find the maximum distance observed in this run
          double maxDistanceThisRun = lcc.stream()
                                          .filter(v -> v.distance != Double.MAX_VALUE)
                                          .mapToDouble(v -> v.distance)
                                          .max()
                                          .orElse(0);
          // Update LmaxDijkstra if the maximum distance found in this run exceeds the current LmaxDijkstra
          LmaxDijkstra = Math.max(LmaxDijkstra, maxDistanceThisRun);
      }

      // Max degree can be calculated once after all runs since it doesn't depend on the start vertex
      maxDegreeInLSP_Dijkstra = lcc.stream()
                                   .mapToInt(v -> v.neighbors.size())
                                   .max()
                                   .orElse(0);
      
 	}
  
    private double getDistance(Vertex u, Vertex v) {
        return Math.sqrt(Math.pow(u.x - v.x, 2) + Math.pow(u.y - v.y, 2));
    }

    public double AStarFindLMax(){

        List<Vertex> LCC = getLargestComponent();

        Vertex end = null;
        double maxDistance = -1.0;
        for(int i = 0; i < LCC.size(); i++){
            for(int j = i+1; j < LCC.size(); j++){
                if(getDistance(LCC.get(i), LCC.get(j)) > maxDistance){
                    maxDistance = getDistance(LCC.get(i), LCC.get(j));
                    end = LCC.get(j);
                };
            }
        }

        for(Vertex v : LCC){
            v.distanceToDestination = getDistance(v,end);
        }
        List<Vertex> S = new ArrayList<Vertex>();
        Queue<Vertex> Q = new PriorityQueue<>(new Comparator<Vertex>() {
            @Override
            public int compare(Vertex o1, Vertex o2) {
                if (o1.distanceToDestination+o1.highestCostToNode > o2.distanceToDestination+o2.highestCostToNode){
                    return -1;
                };

                if(o1.distanceToDestination+o1.distanceToDestination < o2.distanceToDestination+o2.highestCostToNode) {
                    return 1;
                }
                return 0;
            }
        });

        for(Vertex v : LCC){
            Q.add(v);
        }

        while(!Q.isEmpty()){
            Vertex u = Q.poll();
            S.add(u);
            for(Vertex v : u.neighbors){
                if(v.highestCostToNode < u.highestCostToNode + getDistance(u, v)){
                    if(v.currentPathToNode.contains(u)) break;
                    v.currentPathToNode = new ArrayList<>(u.currentPathToNode);
                    v.currentPathToNode.add(u);
                    v.highestCostToNode = u.highestCostToNode + getDistance(u, v);
                    if(S.contains(v)){
                        S.remove(v);
                        Q.add(v);
                    }else{
                        Q.remove(v);
                        Q.add(v);
                    }
                }
            }
        }

        return end.currentPathToNode.size();
    }

    public double BeamSearchFindLMax(int initialBeamWidth, double beamWidthIncrement) {
        List<Vertex> LCC = getLargestComponent();

        // Find the end vertex (destination)
        Vertex end = null;
        double maxDistance = -1.0;
        for (int i = 0; i < LCC.size(); i++) {
            for (int j = i + 1; j < LCC.size(); j++) {
                double distance = getDistance(LCC.get(i), LCC.get(j));
                if (distance > maxDistance) {
                    maxDistance = distance;
                    end = LCC.get(j);
                }
            }
        }

        // Initialize beam search with vertices with highest distance to destination
        int beamWidth = initialBeamWidth;
        Set<Vertex> visited = new HashSet<>();
        Queue<Vertex> queue = new PriorityQueue<>(Comparator.comparingDouble((Vertex o) ->
                o.highestCostToNode + o.distanceToDestination).reversed());

        for (int i = 0; i < beamWidth && i < LCC.size(); i++) {
            Vertex startVertex = LCC.get(i);
            startVertex.highestCostToNode = 0.0; // Start vertex has cost 0
            startVertex.currentPathToNode = new ArrayList<>(); // Initialize path
            startVertex.distanceToDestination = getDistance(startVertex, end); // Heuristic
            queue.add(startVertex);
        }

        double maxPathLength = 0.0;
        while (!queue.isEmpty()) {
            Vertex currentVertex = queue.poll();
            visited.add(currentVertex);

            if (currentVertex.currentPathToNode.size() > maxPathLength) {
                maxPathLength = currentVertex.currentPathToNode.size();
            }

            // Explore neighbors of current vertex
            for (Vertex neighbor : currentVertex.neighbors) {
                if (!visited.contains(neighbor) && !currentVertex.currentPathToNode.contains(neighbor)) {
                    neighbor.currentPathToNode = new ArrayList<>(currentVertex.currentPathToNode);
                    neighbor.currentPathToNode.add(currentVertex);
                    neighbor.highestCostToNode = currentVertex.highestCostToNode + 1; // Increment cost
                    neighbor.distanceToDestination = getDistance(neighbor, end); // Update heuristic
                    queue.add(neighbor);
                }
            }

            // Increase beam width dynamically
            if (queue.isEmpty()) {
                beamWidth += beamWidthIncrement;
                for (int i = 0; i < beamWidth && i < LCC.size(); i++) {
                    Vertex startVertex = LCC.get(i);
                    if (!visited.contains(startVertex)) {
                        startVertex.highestCostToNode = 0.0; // Start vertex has cost 0
                        startVertex.currentPathToNode = new ArrayList<>(); // Initialize path
                        startVertex.distanceToDestination = getDistance(startVertex, end); // Heuristic
                        queue.add(startVertex);
                    }
                }
            }
        }

        return maxPathLength;
    }

    public static void main(String[] args) {
        try {
        	//GraphAnalyzer analyzer = new GraphAnalyzer("n50r0.057.edges");
        	//GraphAnalyzer analyzer = new GraphAnalyzer("DSJC500-5.mtx");
            //GraphAnalyzer analyzer = new GraphAnalyzer("inf-euroroad.edges");
            //GraphAnalyzer analyzer = new GraphAnalyzer("inf-power.mtx");
            GraphAnalyzer analyzer = new GraphAnalyzer("n300r0point08VLCC284.edges");
            //GraphAnalyzer analyzer = new GraphAnalyzer("n400r0point064VLCC324.edges");
            //GraphAnalyzer analyzer = new GraphAnalyzer("n500r0point057VLCC400.edges");
            
            analyzer.runAlgorithmsAndCalculateMetrics();
            
            System.out.println("Algorithm\t|VLCC|\tΔ(LCC)\tk(LCC)\tLmax");
            
            // DFS
            // Perform DFS to calculate LmaxDFS and other metrics
            System.out.println("DFS\t\t" + analyzer.calculateVLCC() + "\t" + analyzer.maxDegreeInLSP_DFS + "\t" + 
                               analyzer.calculateAverageDegreeLCC() + "\t" + analyzer.LmaxDFS);
            
            // Dijkstra
            // Perform Dijkstra to calculate LmaxDijkstra and other metrics
            System.out.println("Dijkstra\t" + analyzer.calculateVLCC() + "\t" + analyzer.maxDegreeInLSP_Dijkstra + "\t" + 
                               analyzer.calculateAverageDegreeLCC() + "\t" + analyzer.LmaxDijkstra);

            analyzer.LmaxAStar = analyzer.AStarFindLMax();
            // A*
            // Perform A* to calculate LmaxAStar and other metrics
            System.out.println("A*\t\t" + analyzer.calculateVLCC() + "\t" + analyzer.calculateDeltaLCC() + "\t" + 
                               analyzer.calculateAverageDegreeLCC() + "\t" + analyzer.LmaxAStar);

            analyzer.LmaxBeamSearch = analyzer.BeamSearchFindLMax(5,5);
            //Anytime A*
            //Perform AnytimeA* to calculate BeamSearchAStar and other metrics
            System.out.println("BeamSearchAStar\t\t" + analyzer.calculateVLCC() + "\t" + analyzer.calculateDeltaLCC() + "\t" +
                    analyzer.calculateAverageDegreeLCC() + "\t" + analyzer.LmaxBeamSearch);


        } catch (IOException e) {
            System.err.println("Error reading the graph file: " + e.getMessage());
        }
    }
}