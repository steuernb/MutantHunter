package coreClasses;

import java.util.Collections;
import java.util.Comparator;
import java.util.Enumeration;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Vector;

import nlr_parser.MastMotifHitList;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;

import supportClasses.BioSequence;
import supportClasses.MPileupLine;

public class NLR_Contig {

	String contigName;
	Hashtable<String,Hashtable<Integer,MPileupLine>> snps;
	Hashtable<String,int[]> coverages;
	
	/**
	 * regions that should be included in analysis. The key is a pair of int denoting start and end. The integer is the reading frame, which should be -1,-2,-3, 1,2,3. If it is 0 it cannot be used.
	 */
	HashSet<int[]> regions;  
	
	
	String wt_sequence;
	MastMotifHitList motif_list_fw; 
	MastMotifHitList motif_list_rv;
	
	
	public final String wildtype = "wildtype";
	String report;
	
	
	
	
	
	
	
	
	
	
	
	/* ***************************************************
	 * ******   Constructors        **********************
	 * ***************************************************
	 */
	
	
	public NLR_Contig (String contigName){
		this.contigName = contigName;
		this.snps = new Hashtable<String,Hashtable<Integer,MPileupLine>>();
		this.coverages = new Hashtable<String,int[]>();
		this.regions = new HashSet<int[]>();
	}
	
	
	public NLR_Contig (Element xmlElement){
		this.contigName = xmlElement.getAttribute("name");
		this.snps = new Hashtable<String,Hashtable<Integer, MPileupLine>>();
		this.coverages = new Hashtable<String,int[]>();
		this.regions = new HashSet<int[]>();
		
		
		
		//SNPs
		Element snpsElement = (Element) xmlElement.getElementsByTagName("SNPs").item(0);
		if( snpsElement != null){
			NodeList nodeListSNPs = snpsElement.getElementsByTagName("SNP");
			for( int i = 0; i< nodeListSNPs.getLength(); i++){
				Element snpElement = (Element) nodeListSNPs.item(i);
				String mutantLine = snpElement.getAttribute("name");
				int position = Integer.parseInt(snpElement.getAttribute("position"));
				MPileupLine mpileupline = new MPileupLine(snpElement.getAttribute("value").replaceAll("&#9;","\t")  );
				
				Hashtable<Integer, MPileupLine> h = new Hashtable<Integer, MPileupLine>();
				if(this.snps.containsKey(mutantLine)){
					h = this.snps.get(mutantLine);
				}
				h.put(new Integer(position), mpileupline);
				this.snps.put(mutantLine, h);
				
				
			}
		}
		
		//Regions
		Element regionsElement = (Element) xmlElement.getElementsByTagName("Regions").item(0);
		if(regionsElement != null){
			NodeList nodeList = regionsElement.getElementsByTagName("Region");
			for( int i = 0; i < nodeList.getLength(); i++){
				Element regionElement = (Element) nodeList.item(i);
				int[] region = {Integer.parseInt(regionElement.getAttribute("start")), Integer.parseInt(regionElement.getAttribute("end"))  };
				
				if(regionElement.hasAttribute("frame")){
					int[] region2 = {region[0], region[1],Integer.parseInt(regionElement.getAttribute("frame")) };
					region = region2;
				}	
				this.regions.add(region);
			}
			
		}
		
		
		
		//coverages
		Element coveragesElement = (Element) xmlElement.getElementsByTagName("Coverages").item(0);
		if(coveragesElement != null){
			NodeList nodeList = coveragesElement.getElementsByTagName("Coverage");
			for( int i = 0; i< nodeList.getLength(); i++){
				Element coverageElement = (Element) nodeList.item(i);
				String mutantLine = coverageElement.getAttribute("name");
				String[] a = coverageElement.getAttribute("value").split(",");
				int[] cov = new int[ a.length];
				for( int j = 0; j< a.length; j++){
					
					cov[j] = Integer.parseInt(a[j]);
				}
				this.coverages.put(mutantLine, cov);
			}
			
		}
		
		
		//sequence
		Element sequenceElement = (Element) xmlElement.getElementsByTagName("sequence").item(0);
		if( sequenceElement != null){
			this.wt_sequence = sequenceElement.getAttribute("wt_dna");
		}
		
		
		//MastMotifHitLists
		Element motifListsElement = (Element) xmlElement.getElementsByTagName("MastMotifHitLists").item(0);
		if( motifListsElement != null){
			NodeList nodeList = motifListsElement.getElementsByTagName("MastMotifHits");
			for( int i = 0; i< nodeList.getLength(); i++){
				Element listElement = (Element) nodeList.item(i);
				if( !listElement.hasAttribute("strand")){
					continue;
				}
				if(listElement.getAttribute("strand").equalsIgnoreCase("fw")){
					this.motif_list_fw = new MastMotifHitList(listElement);
				}else if(listElement.getAttribute("strand").equalsIgnoreCase("rv")){
					this.motif_list_rv = new MastMotifHitList(listElement);
				}
			}
		}
	}
	
	
	
	
	
	
	/* ********************************************* *
	 * ********ADD and SET Methods      ************ * 
	 * ********************************************* */
	
	
	
	
	public void addRegion(int[] a){
			this.regions.add(a);
	}
	

	/**
	 * delete redundancy in regions
	 */
	public void consolidateRegions(){
		int maxIndex = 0;
		if(this.wt_sequence != null){
			maxIndex = this.wt_sequence.length();
		}else{
			for( Iterator<int[]> iterator = this.regions.iterator(); iterator.hasNext();){
				int[] region = iterator.next();
				if(region[1]>maxIndex){
					maxIndex = region[1];
				}
			}
		}
		
		int[][] a = new int[7][maxIndex];
		
		for( Iterator<int[]> iterator = this.regions.iterator(); iterator.hasNext();){
			int[] region = iterator.next();
			int frame = 0;
			if( region.length>2){
				frame = region[2];
			}
			int j =0;
			if( frame ==-3){j=1;}
			if( frame ==-2){j=2;}
			if( frame ==-1){j=3;}
			if( frame ==1){j=4;}
			if( frame ==2){j=5;}
			if( frame ==3){j=6;}
			for( int i = region[0]-1; i< region[1]; i++ ){
				
				a[j][i]++;
			}
		}
		this.regions.clear();
		int start = -1;
		for( int j = 0; j< 7;j++){
			int frame = 0;
			if( j==1){frame = -3;}
			if( j==2){frame = -2;}
			if( j==3){frame = -1;}
			if( j==4){frame = 1;}
			if( j==5){frame = 2;}
			if( j==6){frame = 3;}
			
			for(int i = 0; i< maxIndex; i++){
				if( a[j][i]>0){
					if( start == -1){
						start = i+1;
					}
				}else{
					if( start != -1){
						int[] region = {start, i+1, frame};
						this.regions.add(region);
						start = -1;
					}
				}
				
			}
		}
		
		
		
	}
	
	
	public void addSNP(String mutantLine, MPileupLine line){
	
		
		
		int pos = line.getPosition();
		
		if( this.snps.containsKey(mutantLine)){
			this.snps.get(mutantLine).put(new Integer(pos),line);
		}else{
			Hashtable<Integer, MPileupLine> h = new Hashtable<Integer, MPileupLine>();
			h.put(new Integer(pos), line);
			this.snps.put(mutantLine, h);
		}
	}
	
	public void addWildtypeSNP(MPileupLine line){
		addSNP(wildtype, line);
	}
	
	
	
	public void addCoverage(String mutantLine, int[] coverage){
		
		this.coverages.put(mutantLine, coverage);
	}
	
	public void addWiledtypeCoverage(int[] coverage){
		addCoverage(wildtype, coverage);
	}
	
	
	
	/**
	 * 
	 * Merge another NLR_Contig into this one. Only works if both contigs have the same name and if both have a sequence that sequence has to be identical.
	 * if the new contig has elements that are not present in the current, they will be added to the current.
	 * 
	 * 
	 * @param contig2
	 * 		The NLR_Contig to merge into this one.
	 * @return
	 * 		true if the merging worked smoothly.
	 * 		false if the contig names are inconsistent or both have sequences and the sequences do not match.
	 */
	public boolean mergeContig(NLR_Contig contig2){
		if( this.contigName.equalsIgnoreCase(contig2.contigName)){
			if( this.wt_sequence != null && contig2.getSequence() != null){
				if( !this.wt_sequence.equalsIgnoreCase(contig2.getSequence())){
					return false;
				}
			}
			if( this.wt_sequence == null ){
				this.wt_sequence = contig2.getSequence();
			}
			if( this.motif_list_fw == null){
				this.motif_list_fw = contig2.motif_list_fw;
			}
			if( this.motif_list_rv == null){
				this.motif_list_rv = contig2.motif_list_rv;
			}
			
			
			for(Iterator<int[]> iterator = contig2.regions.iterator(); iterator.hasNext();){
				int[] key = iterator.next();
				this.addRegion(key);
			}
			this.consolidateRegions();
			
			for(Enumeration<String> myenum = contig2.coverages.keys(); myenum.hasMoreElements();){
				String mutantLine = myenum.nextElement();
				if( !this.coverages.containsKey(mutantLine)){
					this.coverages.put(mutantLine, contig2.coverages.get(mutantLine));
				}
			}
			
			for(Enumeration<String> myenum = contig2.snps.keys(); myenum.hasMoreElements();){
				String mutantLine = myenum.nextElement();
				if(!this.snps.containsKey(mutantLine)){
					this.snps.put(mutantLine, contig2.snps.get(mutantLine));
				}else{
					Hashtable<Integer, MPileupLine> h = this.snps.get(mutantLine);
					for(Enumeration<Integer> myenum2 = contig2.snps.get(mutantLine).keys(); myenum2.hasMoreElements();){
						Integer i = myenum2.nextElement();
						if(!h.containsKey(i)){
							h.put(i, contig2.snps.get(mutantLine).get(i));
						}
					}
					this.snps.put(mutantLine, h);
				}
				
			}
			
			
			
			return true;
		}
			
		
		return false;
	}
	
	
	
	
	
	
	
	
	
	
	
	
	public void setSequence(String sequence){
		this.wt_sequence = sequence;
	}
	
	public void setMotifsFw(MastMotifHitList list){
		this.motif_list_fw = list;
	}
	public void setMotifsRv(MastMotifHitList list){
		this.motif_list_rv = list;
	}
	
	public void setRegions(HashSet<int[]> regions){
		this.regions = regions;
	}
	
	
	
	
	
	
	
	/* ********************************************* *
	 * ********        GET Methods      ************ * 
	 * ********************************************* */
	
	
	
	
	public String getContigName(){
		return this.contigName;
	}
	
	
	
	
	public String getReport(){
		return this.report;
	}
	
	public String getSequence(){
		return this.wt_sequence;
	}
	
	public int[] getCoverage(String mutantLine){
		return this.coverages.get(mutantLine);
	}
	
	
	public Hashtable<String,Hashtable<Integer,MPileupLine>> getSnps(){
		return this.snps;
	}
	
	public HashSet<int[]> getRegions(){
		return this.regions;
	}
	public String getRegionReport(){
		//System.err.println(this.getContigName() + "\t" + this.regions.size());
		String s = "";
		Vector<int[]> v = new Vector<int[]>();
		for(Iterator<int[]> iterator = this.regions.iterator(); iterator.hasNext(); ){
			v.add(iterator.next());
		}
		
		Collections.sort(v, new Comparator<int[]>(){
								public int compare(int[] a1, int[]a2){
									if(a1[0] < a2[0]){
										return -1;
									}
									if(a1[0] > a2[0]){
										return 1;
									}
									if(a1[0] == a2[0]){
										if(a1[1] < a2[1]){
											return -1;
										}
										if(a1[1] > a2[1]){
											return 1;
										}
									}
									
									return 0;
								} 
							}     
						);
		
		for(Iterator<int[]> iterator = v.iterator(); iterator.hasNext(); ){
			int[] a = iterator.next();
			s = s+"; "+a[0] + "-" + a[1];
			if(a.length>2){
				s = s + ","+a[2];
			}
		}
		if(s.length()>0){
			s = s.substring(1);
		}
		
		return s;
	}
	
	public Vector<String> getMutantLineNames(){
		
		//sort the lines
		HashSet<String> hashset = new HashSet<String>();
		for(Enumeration<String> myenum = snps.keys(); myenum.hasMoreElements();){
			hashset.add(myenum.nextElement());
		}
		for(Enumeration<String> myenum = coverages.keys(); myenum.hasMoreElements();){
			hashset.add(myenum.nextElement());
		}
		hashset.remove(wildtype);
		Vector<String> mutantLines = new Vector<String>();
		for(Iterator<String> iterator = hashset.iterator();iterator.hasNext();){
			mutantLines.add(iterator.next());
		}
		Collections.sort(mutantLines);
		
		return mutantLines;
	}
	
	
	/**
	 * Gives the length of the contig. If the sequence is added to the object, the return value is the length of that sequence.
	 * Otherwise it tryes to get a length estimation from the coverages.
	 * 
	 * @return
	 */
	public int getLength(){
		if(this.wt_sequence != null){
			return this.wt_sequence.length();
		}
		
		int length = 0;
	
		for(Enumeration<int[]> myenum = this.coverages.elements();myenum.hasMoreElements();){
			int l = myenum.nextElement().length;
			
			if( l> length){
				length = l;
			}
		}
		
		return length;
	}
	
	public int getMedianCoverage(String mutantLine){
		if( this.coverages.get(mutantLine) ==null){
			return 0;
		}
		Vector<Integer> v = new Vector<Integer>();
		for( int i = 0; i< this.coverages.get(mutantLine).length; i++){
			if( this.isInRegion(i+1)){
				v.add(this.coverages.get(mutantLine)[i]);
			}
		}
		Collections.sort(v);
		if( v.size()>0){
			return v.get(v.size()/2).intValue();
		}else{
			return 0;
		}
		
	}
	
	public int getWildtypeCoverage(){
		return this.getMedianCoverage(MutantHunter.wildtype);
	}
	
	
	
	
	/* ********************************************* *
	 * ********        Export           ************ * 
	 * ********************************************* */
	
	

	
	
	public Element getXMLElement(Document dom){
		Element contigElement = dom.createElement("NBLRR_Contig");
		
		
		
		
		contigElement.setAttribute("name",this.getContigName());
		contigElement.setAttribute("length", this.getLength()+"");
		
		
		Element snpsElement = dom.createElement("SNPs");
		contigElement.appendChild(snpsElement);
		for(Enumeration<String> myenum2 = this.snps.keys(); myenum2.hasMoreElements();){
			String mutantLine = myenum2.nextElement();
			Hashtable<Integer, MPileupLine> snps = this.snps.get(mutantLine);
			for(Enumeration<Integer> myenum3 = snps.keys(); myenum3.hasMoreElements();){
				Integer position = myenum3.nextElement();
				
				Element snpElement = dom.createElement("SNP");
				snpElement.setAttribute("name", mutantLine);
				snpElement.setAttribute("position", position.intValue()+"");
				
				String mpileupLine = snps.get(position).toString();
				mpileupLine = mpileupLine.replaceAll(">", ";"); //avoid these symbols in XML by changing quality value by 1. 
				mpileupLine = mpileupLine.replaceAll("<", ";");
				mpileupLine = mpileupLine.replaceAll("\"", "!");
				snpElement.setAttribute("value", mpileupLine);
				snpsElement.appendChild(snpElement);
			}
		}
		
		
		
		if(this.regions!= null){
			Element regionsElement = dom.createElement("Regions");
			contigElement.appendChild(regionsElement);
			for(Iterator<int[]> iterator = this.regions.iterator(); iterator.hasNext();){
				int[] region = iterator.next();
				Element regionElement = dom.createElement("Region");
				regionElement.setAttribute("start", region[0]+"");
				regionElement.setAttribute("end", region[1]+"");
				if( region.length>2){
					regionElement.setAttribute("frame", region[2]+"");
				}	
				regionsElement.appendChild(regionElement);
			}
		}
		
		
		Element coveragesElement = dom.createElement("Coverages");
		contigElement.appendChild(coveragesElement);
		for( Enumeration<String> myenum2 = this.coverages.keys(); myenum2.hasMoreElements();){
			String mutantLine = myenum2.nextElement();
			int[] cov = this.coverages.get(mutantLine);
			String s = "";
			for( int i= 0; i< cov.length; i++){
				s = s +"," +cov[i];
			}
			Element coverageElement = dom.createElement("Coverage");
			coverageElement.setAttribute("name", mutantLine);
			coverageElement.setAttribute("value", s.substring(1));
			coveragesElement.appendChild(coverageElement);
		}
		
		
		Element sequenceElement = dom.createElement("sequence");
		sequenceElement.setAttribute("wt_dna", this.wt_sequence);
		contigElement.appendChild(sequenceElement);
		
		 	
		Element motifListsElement = dom.createElement("MastMotifHitLists");
		contigElement.appendChild(motifListsElement);
		if( this.motif_list_fw != null){
			MastMotifHitList list = this.motif_list_fw;
			Element listElement = list.getXMLElement(this.contigName+"+", dom);
			listElement.setAttribute("strand", "fw");
			motifListsElement.appendChild(listElement);
		}
		if( this.motif_list_rv != null){
			MastMotifHitList list = this.motif_list_rv;
			Element listElement = list.getXMLElement(this.contigName+"-", dom);
			listElement.setAttribute("strand", "rv");
			motifListsElement.appendChild(listElement);
		}	
	
		
		
		return contigElement;
	}
	
	
	
	
	
	
	
	
	
	
	
	/* ********************************************* *
	 * ********           FILTERS       ************ * 
	 * ********************************************* */
	
	
	
	
	
	
	public void filterCommonSNPs(Vector<String> mutantLines){
		HashSet<Integer> usedSNPs = new HashSet<Integer>();
		HashSet<Integer> removable = new HashSet<Integer>();
		
		for( Enumeration<String> myenum1 = mutantLines.elements(); myenum1.hasMoreElements();){
			String mutantLine = myenum1.nextElement();
			if( !this.snps.containsKey(mutantLine)){
				continue;
			}
			for(Enumeration<Integer> myenum2 =  snps.get(mutantLine).keys(); myenum2.hasMoreElements();){
				Integer integer = myenum2.nextElement();
				if(usedSNPs.contains(integer)){
					removable.add(integer);
				}
				usedSNPs.add(integer);
			}
		}
		
		
		for(Iterator<Integer> iterator = removable.iterator(); iterator.hasNext();){
			Integer integer = iterator.next();
			for( Enumeration<String> myenum = mutantLines.elements(); myenum.hasMoreElements();){
				String mutantLine = myenum.nextElement();
				try{
					snps.get(mutantLine).remove(integer);
					//System.out.println("SNP removed: "+mutantLine + "\t" + integer.intValue());
				}catch (NullPointerException e){}	
			}
		}
	}
	
	/**
	 * removes a SNP from a mutant line if the same SNP is also present in the wildtype and is thus rather an assembly mistake or a heterozygous position than an EMS mutation. 
	 * 
	 * If there is a huge difference in the allele frequency in the wildtype and the allele frequency in the mutantLine this could also be an EMS mutation combined with an assembly mistake. 
	 * This can be detected by applying the maxDifferenceInAlleleFrequency threshold.
	 * 
	 * @param mutantLines
	 * @param maxDifferenceInAlleleFrequency
	 */
	public void filterWildtypeSNPs(Vector<String> mutantLines, double maxDifferenceInAlleleFrequency){
		if( snps.get(wildtype) != null){
			for(Enumeration<Integer> myenum1 = snps.get(this.wildtype).keys(); myenum1.hasMoreElements(); ){
				Integer integer = myenum1.nextElement();
				MPileupLine wt_snp = snps.get(this.wildtype).get(integer);
				for(Enumeration<String> myenum2 = snps.keys(); myenum2.hasMoreElements();){
					String mutantLine = myenum2.nextElement();
					if(mutantLine.equalsIgnoreCase(this.wildtype)){
						continue;
					}else{
						if( this.snps.get(mutantLine).containsKey(integer)){
							MPileupLine snp = snps.get(mutantLine).get(integer);
							if( Math.abs(snp.getReferenceAlleleFrequency() - wt_snp.getReferenceAlleleFrequency() ) < maxDifferenceInAlleleFrequency ){
								snps.get(mutantLine).remove(integer);
								//System.out.println("Removed WT comparison: " + mutantLine +"\t" + this.getContigName() + "\t" + integer.intValue() +"\t" + snp.getReferenceAlleleFrequency());
							}
							
						}
						
					}
				}
			}
		}
		
	}
	
	
	
	
	
	/**
	 * delete every SNP that is outside defined regions. If no regions are set, i.e. field regions is either null or empty, no SNP is deleted.
	 * 
	 */
	public void filterForRegions(){
		if(this.regions!=null && regions.size()>0){
			Hashtable<String,Hashtable<Integer, MPileupLine>> newSNPs = new Hashtable<String,Hashtable<Integer, MPileupLine>>();
			
			for(Enumeration<String> myenum1 = this.snps.keys(); myenum1.hasMoreElements();){
				String mutantLine = myenum1.nextElement();
				Hashtable<Integer, MPileupLine> h= new Hashtable<Integer, MPileupLine>();
				
				for(Enumeration<Integer> myenum2 = this.snps.get(mutantLine).keys(); myenum2.hasMoreElements();){
					int position = myenum2.nextElement().intValue();
					
					for(Iterator<int[]> iterator = this.regions.iterator(); iterator.hasNext();){
						int[] region = iterator.next();
						if( region[0] <= position && position <= region[1] ){
							h.put(position, this.snps.get(mutantLine).get(position));
						//	System.out.println("SNP kept: "+mutantLine + "\t" + position + "\t" + snps.get(mutantLine).get(position).getReferenceAlleleFrequency());
						}
					}
				}
				
				newSNPs.put(mutantLine, h);
				
			}
			this.snps = newSNPs;
		}
		
	}
	
	
	
	
	
	
	
	
	
	
	
	/* ********************************************* *
	 * ********           ANALYSIS      ************ * 
	 * ********************************************* */
	
	
	public boolean isInRegion(int position){
		if(this.regions== null || this.regions.size() ==0 ){
			return true; 
		}
		
		for(Iterator<int[]> iterator = this.regions.iterator();iterator.hasNext();){
			int[] region = iterator.next();
			if( position >= region[0] && position <= region[1]){
				return true;
			}
		}
		return false;
	}
	
	
	public int countZeroCoveragePositions(String mutantLine, int minWtCoverage){
		
		int[] regionCoverage = new int[wt_sequence.length()];
		for( Iterator<int[]> iterator = this.regions.iterator(); iterator.hasNext();){
			int[] region = iterator.next();
			for( int i = region[0]-1; i< region[1]; i++){
				regionCoverage[i]++;
			}
		}
		
		
		int numZero = 0;
		for( int i = 0; i<regionCoverage.length; i++){
			if( regionCoverage[i] >0 && this.coverages.get(wildtype)!= null && this.coverages.get(wildtype)[i]>minWtCoverage){
				if(this.coverages.get(mutantLine)==null ||this.coverages.get(mutantLine)[i]==0 ){
					numZero++;
				}
				
				
			}
				
			
		}
		
		
		
		return numZero;
			
	}
	
	
	
	
	
	
	
	
	
	public boolean isCandidate(Vector<String> mymutants, 
										int minWildtypeCoverage, 
										double maxReferenceAlleleFrequency, 
										int minCoverageToConsiderSNP,
										int minNumberOfZeroCoveragePositions,
										int minNumberOfTotalMutants, 
										boolean filterSynonymous){
		
		//initialize report
				this.report = this.getContigName() + "\tlength:"+this.getLength() + "\n";
				
				
				//check wildtype median coverage
				int wildTypeCoverage = this.getWildtypeCoverage();
				if(wildTypeCoverage< minWildtypeCoverage ){
					return false;
				}
				
				
				String linereports = "";
				
				
				
				//check snp mutations
				HashSet<String> mutantLinesWithSNPs = new HashSet<String>();
				HashSet<String> mutantLinesWithDeletions = new HashSet<String> ();
				
				
				for(Enumeration<String> myenum1 =mymutants.elements(); myenum1.hasMoreElements();){
					String mutantLine = myenum1.nextElement();
					linereports = linereports + mutantLine + ": ";
					if( mutantLine.equalsIgnoreCase(wildtype)){
						continue;
					}
					
					int zeroCovPos = this.countZeroCoveragePositions(mutantLine, minWildtypeCoverage);
					if( zeroCovPos >= minNumberOfZeroCoveragePositions){
						linereports = linereports + "\tdeletion mutant(" + zeroCovPos + ")";
						mutantLinesWithDeletions.add(mutantLine);
					}
					
					
					if( this.snps.containsKey(mutantLine)){
						String snpreport = new String("");
						for(Enumeration<Integer> myenum2 = this.snps.get(mutantLine).keys(); myenum2.hasMoreElements();){
							Integer key = myenum2.nextElement();
							if(this.snps.get(mutantLine).get(key).getReferenceAlleleFrequency() <= maxReferenceAlleleFrequency && this.snps.get(mutantLine).get(key).getCoverage()>= minCoverageToConsiderSNP){
								
								
								
								//check if nonsynonymous
								Vector<char[]> v = this.getAminoAcids(key, mutantLine);
								boolean nonsynonymous = false;
								
								String aa_report = "";
								if(v.size()>0){
									aa_report = ","+v.get(0)[0]+"->"+v.get(0)[1];
									boolean changed = false;
									for( Iterator<char[]> iterator = v.iterator(); iterator.hasNext();){
										char[] a = iterator.next();
										if( a[0] !=a[1] && !changed){
											changed = true;
											aa_report = ","+a[0]+"->"+a[1];
											nonsynonymous = true;
										}
									}
									if( v.size()>1){
										aa_report = aa_report + "*";
									}
								}
								
								
								if( nonsynonymous  || !filterSynonymous){
									mutantLinesWithSNPs.add(mutantLine);
									
									snpreport = snpreport + ";SNP("+key.intValue()+ "," + snps.get(mutantLine).get(key).getReferenceAlleleFrequency() +
											"," + snps.get(mutantLine).get(key).getRefBase() + "->" + snps.get(mutantLine).get(key).getLargestAlternativeAllele()+aa_report+")";
								}
								
							}
						}
						if( snpreport.length()>0){
							linereports = linereports + "\t" + snpreport.substring(1) ;
						}	
							
					}
					linereports = linereports +"\n";
					
				}
				
				this.report = this.report + "Number Of SNP mutants: " + mutantLinesWithSNPs.size() + "\n";
				this.report = this.report + "Number Of Deletion mutants: " + mutantLinesWithDeletions.size() + "\n";
				this.report = this.report + "Wildtype coverage: "+ wildTypeCoverage + "\n";
				this.report = this.report + "Regions: "+this.getRegionReport()+"\n";
				this.report = this.report + linereports+"\n";
				
				//final result
				if( mutantLinesWithSNPs.size() >= minNumberOfTotalMutants){
					
					return true;
				}
				
				return false;
		
	}
	
	
	
	public boolean isCandidate(int minWildtypeCoverage, 
								double maxReferenceAlleleFrequency, 
								int minCoverageToConsiderSNP,
								int minNumberOfZeroCoveragePositions,
								int minNumberOfTotalMutants,
								boolean filterSynonymous){
		return  isCandidate(this.getMutantLineNames() , minWildtypeCoverage, maxReferenceAlleleFrequency,minCoverageToConsiderSNP,minNumberOfZeroCoveragePositions,minNumberOfTotalMutants,filterSynonymous);
		
	}
	
	
	private Vector<char[]> getAminoAcids(int position, String mutantLine){
		Vector<char[]> v = new Vector<char[]>();
		
		char mutantNucleotide = this.snps.get(mutantLine).get(position).getLargestAlternativeAllele().toCharArray()[0];
		//char wildtypNucleotide = this.snps.get(mutantLine).get(position).getRefBase();
		
		for(Iterator<int[]> iterator = this.regions.iterator(); iterator.hasNext();){
			int[] region = iterator.next();
			if( region.length>2 && region[0]<= position && region[1]>=position){
				boolean forwardDirection = true;
				if(region[2]<0){
					forwardDirection = false;
				}
				int frame = Math.abs(region[2])-1;
				
				
				if(forwardDirection){
					int pos = position -1;
					int positionInTriplet = (pos-frame)%3 ;
					int aa_start = pos - positionInTriplet  ;
					int aa_end   = pos  +3 -positionInTriplet;
					
					String triplet = "NNN";
					try{
						triplet = this.wt_sequence.substring(aa_start, aa_end);
					}catch(StringIndexOutOfBoundsException e){}
					if( positionInTriplet <0){
						triplet = "NNN";
					}
					char wt_aminoacid = new BioSequence("","").translateTriplet(triplet);
					
					char[] tripletChar = triplet.toCharArray();
					try{
						tripletChar[positionInTriplet] = mutantNucleotide;
					}catch(ArrayIndexOutOfBoundsException e){}
					char mutant_aminoacid = new BioSequence("","").translateTriplet(new String(tripletChar));
					
					char[] a = {wt_aminoacid, mutant_aminoacid}; 
					v.add(a);
				}else{
					
					int pos = wt_sequence.length()-position;
					int positionInTriplet = (pos-frame)%3 ;
					int aa_start = pos - positionInTriplet  ;
					int aa_end   = pos  +3 -positionInTriplet;
					
					String triplet = "NNN";
					try{
						triplet = new BioSequence("",this.wt_sequence).getReverseComplementarySequence().substring(aa_start, aa_end);
					}catch(StringIndexOutOfBoundsException e){System.err.println("out of bounds:"+aa_start);}
					if( positionInTriplet <0){
						triplet = "NNN";
					}
					char wt_aminoacid = new BioSequence("","").translateTriplet(triplet);
					
					char[] tripletChar = triplet.toCharArray();
					tripletChar[positionInTriplet] = new BioSequence("", ""+mutantNucleotide).getReverseComplementarySequence().toCharArray()[0];
					
					char mutant_aminoacid = new BioSequence("","").translateTriplet(new String(tripletChar));
					
					char[] a = {wt_aminoacid, mutant_aminoacid}; 
					//System.err.println(new BioSequence("", wt_sequence).getReverseComplementarySequence());
					//System.err.println(position +"\t"+pos+"\t" +region[2]+ "\t"+wildtypNucleotide + "\t" + mutantNucleotide +"\t" + triplet +"\t" +new String(tripletChar)+ "\t" + wt_aminoacid +"\t" + mutant_aminoacid);
					v.add(a);
					
					
					
					
				}
				
				
				
			}
		}
		
		
		
		
		return v;
	}
	
	public String getDeletionReport(Vector<String> mutantLines, int minWtCoverage, int minZeroCoveagePositions){
		String s = "";
		
		int[] covered = new int[wt_sequence.length()];
		for(Iterator<int[]> iterator = this.regions.iterator(); iterator.hasNext();){
			int[] a = iterator.next();
			for( int i = a[0]-1; i<a[1];i++){
				covered[i]++;
			}
		}
		
		int numBasesInRegions = 0;
		
		for( int i = 0; i< covered.length; i++){
			if(covered[i]>0){
				numBasesInRegions++;
			}
		}
		
		
		
		
		for(Enumeration<String> myenum1 = mutantLines.elements(); myenum1.hasMoreElements();){
			String mutantLine = myenum1.nextElement();
			int z = this.countZeroCoveragePositions(mutantLine, minWtCoverage);
			if(numBasesInRegions>100 &&  z >= minZeroCoveagePositions ){
				s = s + mutantLine +"\t" + this.getContigName() + "\tdeletionmutant\t"+numBasesInRegions + "\t"+ z + "\n";
			}
			
			
		}	
		
		return s;
	}
	
	public String getSNPReport(Vector<String> mutantLines, int minWtCoverage, int minMtCoverage){
		String s = "";
		
		for(Enumeration<String> myenum1 = mutantLines.elements(); myenum1.hasMoreElements();){
			String mutantLine = myenum1.nextElement();
			if( this.snps.get(mutantLine)!= null){
				for(Enumeration<Integer>myenum2 = this.snps.get(mutantLine).keys(); myenum2.hasMoreElements();){
					int position = myenum2.nextElement();
					MPileupLine line = this.snps.get(mutantLine).get(position);
					double raf = line.getReferenceAlleleFrequency();
					
					int wtCoverage = 0;
					if(this.coverages.containsKey(wildtype)){
						wtCoverage = this.coverages.get(this.wildtype)[position-1];
					}
							
					if(raf<0.1){
						Vector<char[]> v = this.getAminoAcids(position, mutantLine);
						boolean nonsynonymous = false;
						char wta ='X';
						char mta ='X';
						for(Enumeration<char[]> myenum3 = v.elements(); myenum3.hasMoreElements();){
							char[] c = myenum3.nextElement();
							wta=c[0];
							mta=c[1];
							if(c[0]!= c[1]){
								nonsynonymous = true;
								break;
							}
						}
						if( wtCoverage>=minWtCoverage && line.getCoverage()>= minMtCoverage){
							s = s + mutantLine +"\t"+line.getChromosome()+ "\t" + position +"\t"+ line.getRefBase() + "\t" + line.getLargestAlternativeAllele() +"\t"+wtCoverage+"\t"+line.getCoverage()+"\t" + raf + "\t" +wta +"\t" +mta+"\t"+ nonsynonymous +"\n";
						}
						
					}	
				}
			}	
		}
		
		
		
		
		
		return s;
	}
	
	
	/**
	 * @deprecated
	 * 
	 * @param minWildtypeCoverage
	 * @param maxReferenceAlleleFrequency
	 * @param minCoverageToConsiderSNP
	 * @param minNumberOfZeroCoveragePositions
	 * @param minNumberOfTotalMutants
	 * @return
	 */
	public boolean isCandidate_old(	int minWildtypeCoverage, 
								double maxReferenceAlleleFrequency, 
								int minCoverageToConsiderSNP,
								int minNumberOfZeroCoveragePositions,
								int minNumberOfTotalMutants){
		
		this.report = new String(this.contigName);
		String linereports = new String("");
		
		if(this.wt_sequence != null){
			this.report = this.report + "\t" + this.wt_sequence.length();
		}
		
		System.out.print(this.contigName + "\t" );
		
		HashSet<String> numberOfSNPMutants = new HashSet<String>(); //this remembers the lines with at least one SNP in the contig
		HashSet<String> numberOfDeletionMutants = new HashSet<String>(); //this remembers the lines with a deletion in contig
		for(Enumeration<String> myenum1 = this.getMutantLineNames().elements(); myenum1.hasMoreElements();){
			String mutantLine = myenum1.nextElement();
			
			linereports = linereports + mutantLine+": ";
			
			int numZeroCoveragePositions = this.countZeroCoveragePositions(mutantLine, minWildtypeCoverage);
			if(numZeroCoveragePositions >= minNumberOfZeroCoveragePositions){
				numberOfDeletionMutants.add(mutantLine);
				linereports = linereports + "deletion mutant(" + numZeroCoveragePositions + ")";
			}
			
			
			if( this.snps.containsKey(mutantLine)){
				String snpreport = new String("");
				for(Enumeration<Integer> myenum2 = this.snps.get(mutantLine).keys(); myenum2.hasMoreElements();){
					Integer key = myenum2.nextElement();
					if(this.snps.get(mutantLine).get(key).getReferenceAlleleFrequency() <= maxReferenceAlleleFrequency){
						numberOfSNPMutants.add(mutantLine);
						snpreport = snpreport + ";SNP("+key.intValue()+ "," + snps.get(mutantLine).get(key).getReferenceAlleleFrequency() +
								"," + snps.get(mutantLine).get(key).getRefBase() + "->" + snps.get(mutantLine).get(key).getLargestAlternativeAllele();
					}
				}
				linereports = linereports + snpreport.substring(1) +"\n";
			}
			
			
		}
		System.out.println(numberOfSNPMutants.size() + "\t" + numberOfDeletionMutants.size() + "\t" );
		this.report = this.report + "\n" + "Number of SNP mutant lines: " + numberOfSNPMutants.size() + "\n"
								+ "Number of deletion mutant lines: " + numberOfDeletionMutants.size() + "\n" + linereports;
		
		if(numberOfSNPMutants.size() >= minNumberOfTotalMutants ){
			return true;
		}
		
		
		return false;
	}
	
	
	
	
	
	
}
