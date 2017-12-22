import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

public class Enrichment {

	public String outputfile = "output/pvals_test_full.txt";
	public String inLists = "genelist/list_official.tsv";
	public String backgroundFile = "background/archs4humangenes.txt";
	
	public int randomListNumber = 1000;
	public int randomListLength = 300;
	
	public Fisher2 fisher;
	public HashSet<String> backgroundGenes;
	public HashSet<String> allGMTGenes;
	public HashSet<String> allGenes = new HashSet<String>();
	public HashMap<String, HashSet<String>> gmts;
	public HashMap<String, HashSet<String>> geneLists;
	
	public double average = 0;
	public double averageOverlap = 0;
	public int counter = 0;
	public int minSize = 20;
	public int minSizeGMT = 20;
	
	public int threadCount = 8;
	
	public double[][] pvals;
	public String[] keys;
	String[] gmtKey;
	
	public static void main(String[] args) {
		
		long time = System.currentTimeMillis();
		Enrichment enrichment = new Enrichment();
		enrichment.initialize();
		System.out.println("Initializing time: "+(System.currentTimeMillis()-time)/1000+"s");
		
		//enrichment.calculate();
		//enrichment.printPvalue();
		time = System.currentTimeMillis();
		enrichment.doRandom();
		System.out.println("Random runtime: "+(System.currentTimeMillis()-time)/1000+"s");
	}
	
	public Enrichment() {
		fisher = new Fisher2(40000);
	}
	
	public void initialize() {
		backgroundGenes = readGenes("background/archs4humangenes.txt");
		gmts = readGMT("gmts", minSizeGMT);
		geneLists = readGenelists(inLists, minSize);
		System.out.println("Number of gmt lists: "+gmts.size());
		System.out.println("Number of gene lists: "+geneLists.size());
	}
	
	public void calculate() {
		
		ExecutorService executor = Executors.newFixedThreadPool(threadCount);
		
		keys = geneLists.keySet().toArray(new String[0]);
		gmtKey = gmts.keySet().toArray(new String[0]);
		
		long time = System.currentTimeMillis();
		
		pvals = new double[100][gmtKey.length];
		
		for(int i=0; i<100; i++) {
			Runnable worker = new EnrichmentThread(i);
	        executor.execute(worker);
			counter+=gmtKey.length;
		}
		
		executor.shutdown();
        while (!executor.isTerminated()) {}
        System.out.println("Finished all threads");
		
		System.out.println("Elapsed time: "+(System.currentTimeMillis() - time)+"ms");
	}
	
	public void printPvalue() {
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outputfile)));
			
			System.out.println("cols: "+pvals[0].length);
			
			bw.write("listid");
			for(int j=0; j<pvals[0].length; j++) {
				bw.write("\t"+gmtKey[j].replace(",", "")+","+gmts.get(gmtKey[j]).size());
			}
			bw.write("\n");
			
			for(int i=0; i<pvals.length; i++) {
				bw.write(keys[i]+","+geneLists.get(keys[i]).size());
				for(int j=0; j<pvals[0].length; j++) {
					bw.write("\t"+pvals[i][j]);
				}
				bw.write("\n");
			}
			bw.close();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	public class EnrichmentThread implements Runnable {
		
		String[] lid;
		String gid;
		int ip = 0;
		int jp = 0;
		
	    public EnrichmentThread(int _i){
	       lid = geneLists.get(keys[_i]).toArray(new String[0]);
	       ip = _i;
	    }

	    @Override
	    public void run() {
	    	
	    		String[] keys = gmts.keySet().toArray(new String[0]);
	    		
	    		for(int i=0; i<keys.length; i++) {
	    			int overlap = 0;
	    			HashSet<String> temp = gmts.get(keys[i]);
	    			
	    			for(int j=0; j<lid.length; j++) {
	    				if(temp.contains(lid[j])){
	        	 			overlap++;
	        	 		}
	    			}
	    			
	    			int numGenelist = lid.length;
	    			int totalBgGenes = backgroundGenes.size();
	    			int totalInputGenes = temp.size();
	    			int numOverlap = overlap;
	    			double oddsRatio = (numOverlap*1.0/(totalInputGenes - numOverlap))/(numGenelist*1.0/(totalBgGenes - numGenelist));
	    			double pvalue = fisher.getRightTailedP(numOverlap,(totalInputGenes - numOverlap), numGenelist, (totalBgGenes - numGenelist));	
	    			pvals[ip][i] = pvalue;
	    			
	    		}
	    }
	}
	
	public HashMap<String, HashSet<String>> readGenelists(String _file, int _minSize){
		
		HashMap<String, HashSet<String>> geneLists = new HashMap<String, HashSet<String>>();
		
		try{
			BufferedReader br = new BufferedReader(new FileReader(new File(_file)));
			String line = "";
			
			//while((line = br.readLine()) != null){
			for(int i=0; i<3000000; i++) {
				
				line = br.readLine();
				
				String[] sp = line.trim().toUpperCase().split("\t");
				String lid = sp[0];
				
				if(geneLists.containsKey(lid)) {
					geneLists.get(lid).add(sp[1]);
				}
				else {
					HashSet<String> genes = new HashSet<String>();
					genes.add(sp[1]);
					geneLists.put(lid, genes);
				}
			}
			br.close();
		}
		catch(Exception e){
			e.printStackTrace();
		}
		
		HashMap<String, HashSet<String>> finalLists = new HashMap<String, HashSet<String>>();
		
		String[] keys = geneLists.keySet().toArray(new String[0]);
		for(String key : keys) {
			geneLists.get(key).retainAll(backgroundGenes);
			
			allGenes.addAll(geneLists.get(key));
			
			if(geneLists.get(key).size() >= minSize) {
				finalLists.put(key, geneLists.get(key));
			}
		}
		
		backgroundGenes.retainAll(allGenes);
		
		return finalLists;
	}
	
	public HashMap<String, HashSet<String>> readGMT(String _path, int _minSize){
		HashMap<String, HashSet<String>> gmt = new HashMap<String, HashSet<String>>();
		allGMTGenes = new HashSet<String>();
		
		File folder = new File(_path);
		File[] listOfFiles = folder.listFiles();
		
		for(File f : listOfFiles){
			try{
				BufferedReader br = new BufferedReader(new FileReader(f));
				String line = "";
				
				while((line = br.readLine()) != null){
					String[] sp = line.trim().toUpperCase().replaceAll(",1.0", "").split("\t");
					String g1 = f.getName()+";"+sp[0];
					
					HashSet<String> genes = new HashSet<String>();
					for(int j=2; j<sp.length; j++){
						genes.add(sp[j]);
					}
					
					genes.retainAll(backgroundGenes);
					
					if(genes.size() >= _minSize) {
						allGMTGenes.addAll(genes);
						gmt.put(g1, genes);
					}
				}
				br.close();
			}
			catch(Exception e){
				e.printStackTrace();
			}
		}
		return gmt;
	}
	
	private HashSet<String> readGenes(String _geneListFile){
		HashSet<String> genes = new HashSet<String>();
		try{
			BufferedReader br = new BufferedReader(new FileReader(new File(_geneListFile)));
			String line = "";
			while((line = br.readLine()) != null){
				genes.add(line.trim().toUpperCase());
			}
			br.close();
		}
		catch(Exception e){
			e.printStackTrace();
		}
		return genes;
	}

	private void doRandom(){

		HashSet<String> genelist = readGenes(backgroundFile);
		genelist.retainAll(allGMTGenes);
		
		String[] genelistarr = genelist.toArray(new String[0]);
		System.out.println("Final genes: "+genelist.size());
		
		Random rn = new Random();
		
		HashMap<String, int[]> rankMap = new HashMap<String, int[]>();
		String[] gmtKey = gmts.keySet().toArray(new String[0]);
		
		double[][] pval = new double[randomListNumber][gmtKey.length];
		int[][] rank = new int[randomListNumber][gmtKey.length];
		double[] avgRank = new double[gmtKey.length];
		double[] std = new double[gmtKey.length];
		
		for(int i=0; i<randomListNumber; i++){
			HashSet<String> randomGenes = new HashSet<String>();
			while(randomGenes.size() < randomListLength){
				randomGenes.add(genelistarr[rn.nextInt(genelistarr.length)]);
			}
			
			String[] rg = randomGenes.toArray(new String[0]);

			for (int j=0; j<gmts.size(); j++) {
				int overlap = 0;
				HashSet<String> gmttemp = gmts.get(gmtKey[j]);
				for(int k=0; k<rg.length; k++) {
					if(gmttemp.contains(rg[k])) {
						overlap++;
					}
				}
				
				int numGenelist = rg.length;
	    			int totalBgGenes = genelist.size();
	    			int totalInputGenes = gmttemp.size();
	    			int numOverlap = overlap;
	    			pval[i][j] = fisher.getRightTailedP(numOverlap,(totalInputGenes - numOverlap), numGenelist, (totalBgGenes - numGenelist));	
	    			
			}
		}
		
		for(int i=0; i<pval.length; i++) {
			rank[i] = getRanksArray(pval[i]);
		}

		try{
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File("output/rank_human.txt")));

			for(int i=0; i<gmtKey.length; i++){
				bw.write(gmtKey[i]);
				for(int j=0; j<rank[0].length; j++){
					bw.write("\t"+rank[j][i]);
				}
				bw.write("\n");
			}

			bw.close();
		}
		catch(Exception e){
			e.printStackTrace();
		}
	}
	
	public static int[] getRanksArray(double[] array) {
	    int[] result = new int[array.length];

	    for (int i = 0; i < array.length; i++) {
	        int count = 0;
	        for (int j = 0; j < array.length; j++) {
	            if (array[j] < array[i]) {
	                count++;
	            }
	        }
	        result[i] = count + 1;
	    }
	    return result;
	}
	
}
