package COMP6651_Project;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

public class GraphAnalyzer {
    static class Vertex {
        int id;

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
    private static final double SCALING_FACTOR = 1000; // Adjust as needed

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
                double distance = Math.sqrt(Math.pow(vertices.get(from).x - vertices.get(to).x, 2) + Math.pow(vertices.get(from).y - vertices.get(to).y, 2));
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
        return allComponents.stream().max(Comparator.comparingInt(List::size)).orElse(new ArrayList<>());
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

    private int calculateDeltaLCC_DFS() {
        List<Vertex> lcc = getLargestComponent(); 
        int maxDegree = 0; 
        LmaxDFS = 0; 

        for (Vertex v : lcc) {
            if (!v.visited) { 
                dfsCalculateDegree(v, maxDegree); 
                maxDegree = Math.max(maxDegree, v.neighbors.size()); // Update after each DFS traversal
            }
        }
        return maxDegree;
    }

    private void dfsCalculateDegree(Vertex v, int pathLength) { 
        v.visited = true;
        v.distance = 0 / SCALING_FACTOR; 

        for (Vertex neighbor : v.neighbors) {
            if (!neighbor.visited) {
                dfsCalculateDegree(neighbor, pathLength + 1); 
            }
        }

        LmaxDFS = Math.max(LmaxDFS, pathLength); // Update after each DFS traversal
    }

    private int calculateDeltaLCC_Dijkstra() {
        List<Vertex> lcc = getLargestComponent(); 
        
        for (Vertex v : lcc) {
            v.distance = Double.MAX_VALUE;
        }
        
        LmaxDijkstra = 0; // Initialize LmaxDijkstra once 
        int overallMaxDegree = 0; // To track the maximum across all traversals
        
        for (Vertex start : lcc) {
            int localMaxDegree = dijkstraCalculateDegree(start); // Get max degree from each run
            overallMaxDegree = Math.max(overallMaxDegree, localMaxDegree); 
        }
        
        return overallMaxDegree; 
    }

    private int dijkstraCalculateDegree(Vertex start) {  // Return maxDegree
        // Standard Dijkstra's implementation
        start.distance = 0;
        PriorityQueue<Vertex> pq = new PriorityQueue<>(Comparator.comparingDouble(v -> v.distance));
        pq.offer(start);
        while (!pq.isEmpty()) {
            Vertex u = pq.poll(); 

            // Note: Mark `u` as visited if needed for your Δ(LCC) calculation
            
            for (Vertex v : u.neighbors) {
                double weight = getDistance(u, v); 
                if (v.distance > u.distance + weight) {
                    v.distance = u.distance + weight;
                    v.predecessor = u;
                    pq.offer(v); 

                    // Update LmaxDijkstra conditionally
                    if (v.distance > LmaxDijkstra) { 
                        LmaxDijkstra = v.distance; 
                    }
                }
            }
        }
        
        int maxDegree = 0; // Initialize maxDegree at the start of each traversal
        int localMaxDegree = 0; 
        for (Vertex neighbor : start.neighbors) {
               localMaxDegree = Math.max(localMaxDegree, neighbor.neighbors.size());  
        }
        return localMaxDegree; 
    }
    
    
    private double getDistance(Vertex u, Vertex v) {
        return u.distances.get(u.neighbors.indexOf(v));
    }

      
    // Δ(LCC) and detailed Lmax calculations require specific implementations based on graph analysis.

    public void generateGeometricGraph(int n, double r){

    }


    public static void main(String[] args) {
        try {
            //GraphAnalyzer analyzer = new GraphAnalyzer("DSJC500-5.mtx");
            //GraphAnalyzer analyzer = new GraphAnalyzer("inf-euroroad.edges");
            //GraphAnalyzer analyzer = new GraphAnalyzer("inf-power.mtx");
            //GraphAnalyzer analyzer = new GraphAnalyzer("n300r0point08VLCC284.edges");
            //GraphAnalyzer analyzer = new GraphAnalyzer("n400r0point064VLCC324.edges");
            GraphAnalyzer analyzer = new GraphAnalyzer("n500r0point057VLCC400.edges");
            System.out.println("Algorithm\t|VLCC|\tΔ(LCC)\tk(LCC)\tLmax");
            // DFS
            // Perform DFS to calculate LmaxDFS and other metrics
            System.out.println("DFS\t\t" + analyzer.calculateVLCC() + "\t" + analyzer.calculateDeltaLCC_DFS() + "\t" + 
                               analyzer.calculateAverageDegreeLCC() + "\t" + analyzer.LmaxDFS);
            // Dijkstra
            // Perform Dijkstra to calculate LmaxDijkstra and other metrics
            System.out.println("Dijkstra\t" + analyzer.calculateVLCC() + "\t" + analyzer.calculateDeltaLCC_Dijkstra() + "\t" + 
                               analyzer.calculateAverageDegreeLCC() + "\t" + analyzer.LmaxDijkstra);
            // A*
            // Perform A* to calculate LmaxAStar and other metrics
            System.out.println("A*\t\t" + analyzer.calculateVLCC() + "\t" + analyzer.calculateDeltaLCC() + "\t" + 
                               analyzer.calculateAverageDegreeLCC() + "\t" + analyzer.LmaxAStar);
        } catch (IOException e) {
            System.err.println("Error reading the graph file: " + e.getMessage());
        }
    }
}
