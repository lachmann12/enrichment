import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.text.DecimalFormat;
import java.util.HashMap;
import java.util.HashSet;

public class Cooccurence {
	
	public String[] genes = new String[0];

	public static void main(String[] args) {
		
		Cooccurence co = new Cooccurence();
		long time = System.currentTimeMillis();
//		HashMap<String, HashSet<Integer>> lists = co.readFile();
		HashMap<Integer, HashSet<String>> lists = co.readFile2();
		System.out.println("Elapsed time reading: "+(System.currentTimeMillis() - time)/1000);
		
		System.out.println("Number of lists: "+lists.size());
		
//		time = System.currentTimeMillis();
//		short[][] overlap = co.getOverlap(lists);
//		System.out.println("Elapsed time overlap: "+(System.currentTimeMillis() - time));
		
		
		time = System.currentTimeMillis();
		short[][] overlap = co.getOverlap2(lists);
		System.out.println("Elapsed time overlap: "+(System.currentTimeMillis() - time)/1000);		
		
		time = System.currentTimeMillis();
		//co.writeTable(overlap, lists);
		co.writeTable2(overlap);
		System.out.println("Elapsed time write: "+(System.currentTimeMillis() - time)/1000);
	}
	
	public HashMap<String, HashSet<Integer>> readFile() {
		HashMap<String, HashSet<Integer>> set = new HashMap<String, HashSet<Integer>>();
		
		try{
			BufferedReader br = new BufferedReader(new FileReader(new File("genelist/list_official.tsv")));
			String line = "";
			
			while((line = br.readLine()) != null){
				String[] sp = line.toUpperCase().trim().split("\t");
				Integer i = Integer.parseInt(sp[0]);
				if (set.containsKey(sp[1])) {
					set.get(sp[1]).add(i);
				}
				else {
					HashSet<Integer> temp = new HashSet<Integer>();
					set.put(sp[1], temp);
				}
			}
			br.close();
			
		}
		catch(Exception e){
			e.printStackTrace();
		}
		
		return set;
	}
	
	public HashMap<Integer, HashSet<String>> readFile2() {
		HashMap<Integer, HashSet<String>> set = new HashMap<Integer, HashSet<String>>();
		
		try{
			BufferedReader br = new BufferedReader(new FileReader(new File("genelist/list_official.tsv")));
			String line = "";
			int counter = 0;
			
			while((line = br.readLine()) != null && counter < 10000000){
				//counter++;
				String[] sp = line.toUpperCase().trim().split("\t");
				Integer i = Integer.parseInt(sp[0]);
				if (set.containsKey(i)) {
					set.get(i).add(sp[1]);
				}
				else {
					HashSet<String> temp = new HashSet<String>();
					temp.add(sp[1]);
					set.put(i, temp);
				}
			}
			br.close();
		}
		catch(Exception e){
			e.printStackTrace();
		}
		return set;
	}
	
	public short[][] getOverlap(HashMap<String, HashSet<Integer>> _set) {
		
		short[][] over = new short[_set.size()][_set.size()];
		
		genes = _set.keySet().toArray(new String[0]);
		
		for(int i=0; i<_set.size(); i++) {
			for(int j=i; j<_set.size(); j++) {
		//for(int i=0; i<10; i++) {
		//	for(int j=i; j<10; j++) {
				HashSet<Integer> temp = new HashSet<Integer>(_set.get(genes[i]));
				temp.retainAll(_set.get(genes[j]));
				over[i][j] = (short)temp.size();
				over[j][i] = (short)temp.size();
			}
		}
		
		return over;
	}
	
	public short[][] getOverlap2(HashMap<Integer, HashSet<String>> _set) {
		
		HashSet<String> allgenes = new HashSet<String>();
		for(Integer i : _set.keySet()) {
			allgenes.addAll(_set.get(i));
		}
		genes = allgenes.toArray(new String[0]);
		System.out.println("Number of genes: "+genes.length);
		short[][] over = new short[genes.length][genes.length];
		System.out.println("matrix created");
		HashMap<String, Integer> gi = new HashMap<String, Integer>();
		
		for(int i=0; i<genes.length; i++) {
			gi.put(genes[i], i);
		}
		
		Integer[] keys = _set.keySet().toArray(new Integer[0]);
		
		for(int i=0; i<keys.length; i++) {
			String[] listGenes = _set.get(keys[i]).toArray(new String[0]);
			for(int j=0; j<listGenes.length; j++) {
				int p1 = gi.get(listGenes[j]);
				for(int k=j; k<listGenes.length; k++) {
					int p2 = gi.get(listGenes[k]);
					over[p1][p2]++;
					if(i!=k) {
						over[p2][p1]++;
					}
				}
			}
		}
		
		return over;
	}
	
	public void writeTable(short[][] _over, HashMap<String, HashSet<Integer>> _set) {

		try{
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File("output/list_official_cooccur.tsv")));
			String[] genes = _set.keySet().toArray(new String[0]);
						
			for(int i=0; i<_over.length; i++) {
				bw.write("\t"+genes[i]);
			}
			bw.write("\n");
			
			for(int i=0; i<_over.length; i++) {
				bw.write(genes[i]);
				for(int j=0; j<_over.length; j++) {
					bw.write("\t"+_over[i][j]);
				}
				bw.write("\n");
			}
			bw.close();
			
			bw = new BufferedWriter(new FileWriter(new File("output/list_official_cooccur_prob.tsv")));
					
			for(int i=0; i<_over.length; i++) {
				bw.write("\t"+genes[i]);
			}
			bw.write("\n");
			
			DecimalFormat df = new DecimalFormat();
			df.setMaximumFractionDigits(3);
			
			for(int i=0; i<_over.length; i++) {
				bw.write(genes[i]);
				for(int j=0; j<_over.length; j++) {
					bw.write("\t"+df.format((_over[i][j]/280000.0)/((_over[i][i]/280000.0)*(_over[j][j]/280000.0))));
				}
				bw.write("\n");
			}
			bw.close();
			
		}
		catch(Exception e){
			e.printStackTrace();
		}
		
	}
	
	public void writeTable2(short[][] _over) {

		try{
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File("output/list_official_cooccur.tsv")));
			
			for(int i=0; i<_over.length; i++) {
				bw.write("\t"+genes[i]);
			}
			bw.write("\n");
			
			for(int i=0; i<_over.length; i++) {
				bw.write(genes[i]);
				for(int j=0; j<_over.length; j++) {
					bw.write("\t"+_over[i][j]);
				}
				bw.write("\n");
			}
			
			bw.close();
			
			bw = new BufferedWriter(new FileWriter(new File("output/list_official_cooccur_prob.tsv")));
					
			for(int i=0; i<_over.length; i++) {
				bw.write("\t"+genes[i]);
			}
			bw.write("\n");
			
			DecimalFormat df = new DecimalFormat();
			df.setMaximumFractionDigits(3);
			
			for(int i=0; i<_over.length; i++) {
				bw.write(genes[i]);
				for(int j=0; j<_over.length; j++) {
					bw.write("\t"+df.format((_over[i][j]/280000.0)/((_over[i][i]/280000.0)*(_over[j][j]/280000.0))));
				}
				bw.write("\n");
			}
			bw.close();
		}
		catch(Exception e){
			e.printStackTrace();
		}
	}

}
