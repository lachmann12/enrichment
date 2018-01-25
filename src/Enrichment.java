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

	public String outputfile = "output/pvals_test_shortlist_10_50.txt";
	public String inLists = "genelist/list_official.tsv";
	public String backgroundFile = "background/archs4humangenes.txt";
	
	public int genelistnumber = 1000;
	public int randomListNumber = 1000;
	public int randomListLength = 300;
	
	public Fisher2 fisher;
	public HashSet<String> backgroundGenes;
	public HashSet<String> allGMTGenes;
	public HashSet<String> allGenes = new HashSet<String>();
	public HashMap<String, HashSet<String>> gmts;
	public HashMap<String, String[]> gmtsarray;
	public HashMap<String, HashSet<String>> geneLists;
	
	public double average = 0;
	public double averageOverlap = 0;
	public int counter = 0;
	
	public int minSize = 10;
	public int maxSize = 50;
	
	public int minSizeGMT = 50;
	
	public int threadCount = 8;
	
	public String[] allgenes;
	
	public double[][] pvals;
	public int[][] overlaps;
	public String[] keys;
	String[] gmtKey;
	
	public static void main(String[] args) {
		
		long time = System.currentTimeMillis();
		Enrichment enrichment = new Enrichment();
		enrichment.initialize();
		System.out.println("Initializing time: "+(System.currentTimeMillis()-time)/1000+"s");
		
		time = System.currentTimeMillis();
		
		enrichment.calculate();
		enrichment.printPvalue();
		
		//enrichment.doRandom();
		System.out.println("Random runtime: "+(System.currentTimeMillis()-time)/1000+"s");
	}
	
	public Enrichment() {
		fisher = new Fisher2(40000);
	}
	
	public void initialize() {
		
		backgroundGenes = readGenes("background/archs4humangenes.txt");
		gmts = readGMT("gmts", minSizeGMT);
		gmtKey = gmts.keySet().toArray(new String[0]);
		gmtsarray = new HashMap<String, String[]>();
		
		for(String gk : gmtKey) {
			gmtsarray.put(gk, gmts.get(gk).toArray(new String[0]));
		}
		
		geneLists = readGenelists(inLists, minSize, maxSize);
		System.out.println("Number of gmt lists: "+gmts.size());
		System.out.println("Number of gene lists: "+geneLists.size());
	}
	
	public void calculate() {
		
		ExecutorService executor = Executors.newFixedThreadPool(threadCount);
		
		keys = geneLists.keySet().toArray(new String[0]);
		gmtKey = gmts.keySet().toArray(new String[0]);
		
		long time = System.currentTimeMillis();
		
		pvals = new double[genelistnumber][gmtKey.length];
		overlaps = new int[genelistnumber][gmtKey.length];
		
		for(int i=0; i<pvals.length; i++) {
			Runnable worker = new EnrichmentThread(i);
	        executor.execute(worker);
			counter += gmtKey.length;
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
	    	
	    		String[] gmtkeys = gmts.keySet().toArray(new String[0]);
	    		
	    		HashSet<String> tempo = geneLists.get(keys[ip]);
	    		
	    		for(int i=0; i<gmtkeys.length; i++) {
	    			int overlap = 0;
	    			
	    			String[] temparr = gmtsarray.get(gmtkeys[i]);
	    			
	    			if(lid.length < temparr.length) {
	    				HashSet<String> temp = gmts.get(gmtkeys[i]);
		    			for(int j=0; j<lid.length; j++) {
		    				if(temp.contains(lid[j])){
		        	 			overlap++;
		        	 		}
		    			}
	    			}
	    			else {
	    				
	    				for(int j=0; j<temparr.length; j++) {
		    				if(tempo.contains(temparr[j])){
		        	 			overlap++;
		        	 		}
		    			}
	    			}
	    			
	    			if(overlap > 0) {
		    			int numGenelist = lid.length;
		    			int totalBgGenes = allGMTGenes.size();
		    			int totalInputGenes = temparr.length;
		    			int numOverlap = overlap;
		    			//double oddsRatio = (numOverlap*1.0/(totalInputGenes - numOverlap))/(numGenelist*1.0/(totalBgGenes - numGenelist));
		    			double pvalue = fisher.getRightTailedP(numOverlap,(totalInputGenes - numOverlap), numGenelist, (totalBgGenes - numGenelist));	
		    			pvals[ip][i] = pvalue;
		    			overlaps[ip][i] = overlap;
	    			}
	    			else {
	    				pvals[ip][i] = 1;
	    			}
	    		}
	    }
	}
	
	public HashMap<String, HashSet<String>> readGenelists(String _file, int _minSize, int _maxSize){
		
		HashMap<String, HashSet<String>> geneLists = new HashMap<String, HashSet<String>>();
		
		try{
			BufferedReader br = new BufferedReader(new FileReader(new File(_file)));
			String line = "";
			
			//while((line = br.readLine()) != null){
			for(int i=0; i<2000000; i++) {
				
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
			
			if(geneLists.get(key).size() >= _minSize &&  geneLists.get(key).size() <= _maxSize) {
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
		
		Random rn = new Random();
		pvals = new double[randomListNumber][gmtKey.length];
		
		HashSet<String> genelist = readGenes(backgroundFile);
		genelist.retainAll(allGMTGenes);
		allgenes = genelist.toArray(new String[0]);
		
		geneLists = new HashMap<String, HashSet<String>>();
		
		System.out.println(allgenes.length);
		
		for(int i=0; i<randomListNumber; i++) {
			HashSet<String> tt = new HashSet<String>();
			while(tt.size() < randomListLength) {
				tt.add(allgenes[rn.nextInt(allgenes.length)]);
			}
			geneLists.put(""+i,  tt);
		}
		
		ExecutorService executor = Executors.newFixedThreadPool(threadCount);
		
		keys = geneLists.keySet().toArray(new String[0]);
		gmtKey = gmts.keySet().toArray(new String[0]);
		
		pvals = new double[randomListNumber][gmtKey.length];
		int[][] rank = new int[randomListNumber][gmtKey.length];
		overlaps = new int[randomListNumber][gmtKey.length];
		
		for(int i=0; i<keys.length; i++) {
			Runnable worker = new EnrichmentThread(i);
	        executor.execute(worker);
		}
		
		executor.shutdown();
		
		double[][] pvalt = new double[pvals[1].length][pvals.length];
		
		for(int i=0; i<pvals.length; i++) {
			for(int j=0; j<pvals[1].length; j++) {
				pvalt[j][i] = pvals[i][j];
			}
		}
		
		for(int i=0; i<pvals.length; i++) {
			rank[i] = getRanksArray(pvals[i]);
		}
		
		System.out.println("ok: "+pvals.length+" - "+pvals[1].length);
		
		try{
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File("output/rank_human.txt")));
			
			for(int i=0; i<gmtKey.length; i++){
				bw.write(gmtKey[i]+"\t"+gmts.get(gmtKey[i]).size());
				for(int j=0; j<rank.length; j++){
					bw.write("\t"+rank[j][i]);
				}
				bw.write("\n");
			}
			bw.close();
		}
		catch(Exception e){
			e.printStackTrace();
		}
		
		try{
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File("output/pval_human.txt")));

			for(int i=0; i<gmtKey.length; i++){
				bw.write(gmtKey[i]+"\t"+gmts.get(gmtKey[i]).size());
				for(int j=0; j<pvals.length; j++){
					bw.write("\t"+pvals[j][i]);
				}
				bw.write("\n");
			}
			bw.close();
		}
		catch(Exception e){
			e.printStackTrace();
		}
		
		try{
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File("output/overlap_human.txt")));

			for(int i=0; i<gmtKey.length; i++){
				bw.write(gmtKey[i]+"\t"+gmts.get(gmtKey[i]).size());
				for(int j=0; j<overlaps.length; j++){
					bw.write("\t"+overlaps[j][i]);
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
