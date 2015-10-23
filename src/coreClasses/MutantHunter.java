package coreClasses;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Collections;
import java.util.Comparator;
import java.util.Enumeration;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.transform.OutputKeys;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerConfigurationException;
import javax.xml.transform.TransformerException;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;


import nlr_parser.MastMotifHitList;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

import supportClasses.BioSequence;
import supportClasses.CLI;
import supportClasses.CLIParseException;
import supportClasses.FastaReader;
import supportClasses.MPileupLine;
import supportClasses.blast.BlastHSP;
import supportClasses.blast.BlastHit;
import supportClasses.blast.BlastIPKReader;
import supportClasses.blast.BlastIteration;
import supportClasses.blast.BlastReader;
import supportClasses.blast.BlastXMLReader;





/**
 * 
 * @version 3.0
 * @author steuernb
 *
 */
public class MutantHunter {

	public static final double version = 3.0;
	public static final String wildtype = "wildtype";
	
	Hashtable<String, NLR_Contig> nlr_contigs;
	//Hashtable<String, Integer> medianCoverages;
	
	HashSet<String> mutantLines;

	
	
	/* ********************************************* *
	 * ********       Constructors      ************ * 
	 * ********************************************* */
	
	
	
	
	
	
	public MutantHunter(  ){
		
		this.nlr_contigs     = new Hashtable<String, NLR_Contig>();
		
		this.mutantLines = new HashSet<String>();
	}
	
	public MutantHunter( File wtXML )throws IOException, ParserConfigurationException, SAXException{
		this.nlr_contigs     = new Hashtable<String, NLR_Contig>();
		//this.medianCoverages = new Hashtable<String, Integer>();
		this.mutantLines = new HashSet<String>();
		this.addXML(wtXML, true);
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	/* ********************************************* *
	 * ********         Getters         ************ * 
	 * ********************************************* */
	public static double getVersion(){
		return version;
	}
	
	public Hashtable<String, NLR_Contig> getNLR_Contigs(){
		return this.nlr_contigs;
	}
	
	/*
	public Hashtable<String,Integer> getMedianCoverages(){
		return this.medianCoverages;
	}
	*/
	
	
	
	
	
	
	
	
	
	
	
	
	/* ********************************************* *
	 * ******** High Level Filters      ************ * 
	 * ********************************************* */
	  
	 /* filter the list of contigs or set regions in contigs*/
	  
	 
	public void setNLRContigList( Hashtable<String,MastMotifHitList> nlrs ){
		
		if(this.nlr_contigs == null){
			this.nlr_contigs = new Hashtable<String, NLR_Contig>();
		}
		
		
		for( Enumeration<String> myenum = nlrs.keys(); myenum.hasMoreElements(); ){
			String key = myenum.nextElement();
			
			String contigName = key.substring(0,key.length()-1);
			
			NLR_Contig contig;
			if( this.nlr_contigs.containsKey(contigName)){
				contig = this.nlr_contigs.get(contigName);
			}else{
				contig = new NLR_Contig(contigName);
			}
			
			String direction = key.substring(key.length()-1);
			
			if(direction.equalsIgnoreCase("+")){
				contig.setMotifsFw(nlrs.get(key));
				if(nlrs.get(key).getDNASequence()!= null){
					contig.setSequence(nlrs.get(key).getDNASequence());
				}
			}
			if(direction.equalsIgnoreCase("-")){
				contig.setMotifsRv(nlrs.get(key));
				
				if( contig.getSequence() == null && nlrs.get(key).getDNASequence() != null){
					contig.setSequence(new BioSequence("", nlrs.get(key).getDNASequence()).getReverseComplementarySequence());
				}
			}
			nlr_contigs.put(contigName, contig);
		}
	}
	
	
	
	public void setContigList(File inputTable)throws IOException{
		BufferedReader in = new BufferedReader(new FileReader(inputTable));
		for (String inputline = in.readLine();inputline != null ;inputline = in.readLine() ) {
			if( inputline.trim().startsWith("#")){
				continue;
			}
			String contigName = inputline.trim().split("\t")[0];
			nlr_contigs.put(contigName, new NLR_Contig(contigName));
		}
		in.close();
	}
	
	public void setContigList(Vector<String> contigNames)throws IOException{
		
		for (Enumeration<String> myenum = contigNames.elements(); myenum.hasMoreElements(); ) {
			
			String contigName =myenum.nextElement();
			nlr_contigs.put(contigName, new NLR_Contig(contigName));
		}
		
	}
	
	/**
	 * 
	 * add the wildtype sequence for every existing nlr_contig. No new contig will be generated.
	 * 
	 * 
	 * @param sequenceFile
	 * 				location of the sequence file. 
	 * @throws IOException
	 */
	public void addSequences(File sequenceFile)throws IOException{
		FastaReader fastaReader = new FastaReader(sequenceFile);
		for (BioSequence seq = fastaReader.readEntry(); seq != null; seq = fastaReader.readEntry()) {
			if( this.nlr_contigs.containsKey(seq.getIdentifier())){
				this.nlr_contigs.get(seq.getIdentifier()).setSequence(seq.getSequence());
			}
		}
		fastaReader.close();
	}
	
	
	
	public void addRegions(File blast, int numHits, boolean addFrame)throws IOException{
		
		boolean isXML = false;
		
		FileInputStream fis = new FileInputStream(blast);
		byte[] bytes = new byte[2];
		fis.read(bytes);
		int head = ((int) bytes[0] & 0xff) | ((bytes[1] << 8) & 0xff00);
		boolean gzip = GZIPInputStream.GZIP_MAGIC == head;
		fis.close();
		
		BufferedReader in;
		if(gzip){
			
			in = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(blast))));
			
			
		}else{
			in = new BufferedReader(new FileReader(blast));
		}
		
		String inputline = in.readLine();
		if(inputline.contains("xml")){
			isXML=true;
		}
		
		in.close();
		
		
		BlastReader reader;
		if(isXML){
			reader = new BlastXMLReader(blast);
		}else{
			reader = new BlastIPKReader(blast);
		}
		
		for(BlastIteration iteration = reader.readIteration(); iteration!= null; iteration= reader.readIteration()){
			int counta = 0;
			
			NLR_Contig contig = this.nlr_contigs.get(iteration.getQueryID());
			if(contig == null){
				continue;
			}
			for(Enumeration<BlastHit> myenum1 = iteration.getHits().elements(); myenum1.hasMoreElements();){
				BlastHit hit = myenum1.nextElement();
				for(Enumeration<BlastHSP> myenum2 = hit.getHSPs().elements(); myenum2.hasMoreElements();){
					BlastHSP hsp = myenum2.nextElement();
					int[] region = {hsp.getQueryStart(), hsp.getQueryEnd(), 0};
					if(addFrame){
						region[2] = hsp.getQueryFrame();
					}
					contig.addRegion(region);
				}
				
				
				counta++;
				if(counta>=numHits){
					break;
				}
			}
			contig.consolidateRegions();
			this.nlr_contigs.put(iteration.getQueryID(), contig);
			
		}
		
		
		reader.close();
		
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	/**
	 * 
	 *  Read a pileup file and add coverage and snps to existing(!) NLR_Contigs. 
	 * 
	 * 
	 * @param mutantLine
	 * 			The name of the line this mapping is from. 	
	 * 	
	 * @param pileup
	 * 			The location of the pileup file.
	 * 
	 * 
	 * @throws IOException
	 */
	public void readPileupFile(String mutantLine, File pileup, int minCoverageForSNP)throws IOException{
		
		
		BufferedReader in;
		
		//check if the fastq is gzipped
		FileInputStream fis = new FileInputStream(pileup);
		byte[] bytes = new byte[2];
		fis.read(bytes);
		int head = ((int) bytes[0] & 0xff) | ((bytes[1] << 8) & 0xff00);
		boolean gzip = GZIPInputStream.GZIP_MAGIC == head;
		fis.close();
		if(gzip){
			in = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(pileup))));
		}else{
			in = new BufferedReader((new FileReader(pileup)));
		}
		
		
		Vector<String> v = new Vector<String>();
		for( String inputline = in.readLine(); inputline != null; inputline = in.readLine()){
			String[] split = inputline.split("\t");
			if( this.nlr_contigs.containsKey(split[0])){
				v.add(inputline);
			}
		}
		in.close();
		//System.out.println("pileup read");
		
		NLR_Contig contig = this.nlr_contigs.get(v.get(0).split("\t")[0]);
		int[] coverage = new int[contig.getLength()];
		for(Enumeration<String> myenum = v.elements(); myenum.hasMoreElements();){
			String inputline = myenum.nextElement();
			String[] split = inputline.split("\t");
			
			if( !split[0].equalsIgnoreCase(contig.getContigName())){
				contig.addCoverage(mutantLine, coverage);
				contig = this.nlr_contigs.get(split[0]);
				coverage = new int[contig.getLength()];
			}
			
			int cov = Integer.parseInt(split[3]);
			int pos = Integer.parseInt(split[1]);
			coverage[pos-1] = cov;
			
			if( cov >= minCoverageForSNP){
				
				Pattern p = Pattern.compile("([atgcATGC]+)");
				Matcher m = p.matcher(split[4]);
				int basecount = 0;
				while(m.find()){
					basecount = basecount + m.group(0).length();
				}
				if(basecount>minCoverageForSNP){
					contig.addSNP(mutantLine, new MPileupLine(inputline));
				}
			}
			
		}
		
		
		
		
		for(Enumeration<String> myenum = this.nlr_contigs.keys(); myenum.hasMoreElements();){
			String key = myenum.nextElement();
			for(Enumeration<String> myenum2 = this.nlr_contigs.get(key).getMutantLineNames().elements(); myenum2.hasMoreElements();){
				mutantLines.add(myenum2.nextElement());
			}
		}
		
		
	}
	
	

	
	
	
	
	
	
	
	
	
	
	

	
	public void findCandidates( File outputFile, int minWildtypeCoverage, 
			double maxReferenceAlleleFrequency,int minCoverageToConsiderSNP, 
			int minNumberOfZeroCoveragePositions,int minNumberOfTotalMutants,
			boolean filterSynonymous)throws IOException{
		
		Vector<String> mutantLines =this.getMutantLines();
		
		findCandidates(mutantLines, outputFile, minWildtypeCoverage, maxReferenceAlleleFrequency, minCoverageToConsiderSNP, minNumberOfZeroCoveragePositions, minNumberOfTotalMutants, filterSynonymous);
		
	}
	
	public void findCandidates(Vector<String> mutantLines, File outputFile, int minWildtypeCoverage, 
								double maxReferenceAlleleFrequency,int minCoverageToConsiderSNP, 
								int minNumberOfZeroCoveragePositions,int minNumberOfTotalMutants,
								boolean filterSynonymous)throws IOException{
		
		
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));

	//	System.err.println("DEBUG: Number of Mutant Lines: "+this.mutantLines.size());
		
		
	//	System.err.println("DEBUG: finding candidates in " + this.getNumberOfContigs() + " contigs");
	//	System.err.println("DEBUG: writing candidates to " + outputFile.getAbsolutePath());
		Vector<String> contigNames = new Vector<String>();
		
		
		for(Enumeration<String> myenum = nlr_contigs.keys(); myenum.hasMoreElements();){
			contigNames.add(myenum.nextElement());
		}
		Collections.sort(contigNames, new Comparator<String>(){ 
													public int compare(String o1, String o2){
														if( !o1.startsWith("contig") || !o2.startsWith("contig")){
															return o1.compareTo(o2);
														}else{
															int i1 = Integer.parseInt(o1.split("_")[1]);
															int i2 = Integer.parseInt(o2.split("_")[1]);
															if(i1<i2){
																return -1;
															}else if(i1 == i2){
																return 0;
															}else{
																return 1;
															}
														}
														
														
													}
										});
		
		
		for(Enumeration<String> myenum = contigNames.elements(); myenum.hasMoreElements();){
			String key = myenum.nextElement();
			NLR_Contig contig = nlr_contigs.get(key);
			//System.out.println(contig.getContigName());
			
			contig.filterForRegions();
			contig.filterWildtypeSNPs(contig.getMutantLineNames(), 0.1);
			contig.filterCommonSNPs(contig.getMutantLineNames());
			if(contig.isCandidate(mutantLines, minWildtypeCoverage, maxReferenceAlleleFrequency, minCoverageToConsiderSNP, minNumberOfZeroCoveragePositions, minNumberOfTotalMutants, filterSynonymous) ){
				out.write(contig.getReport());
				System.out.println(contig.getReport());
				out.newLine();
			}
			
			
		}
		out.close();
		
	}
	
	public void writeDeletionReport(File outputFile, int minWtCoverage, int minZeroCoveragePositions)throws IOException{
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));

		out.write("mutantLine\tcontigName\tdeletionmutant\tregionLength\tdeletionLength\n");
		for(Enumeration<NLR_Contig> myenum = this.nlr_contigs.elements(); myenum.hasMoreElements();){
			NLR_Contig contig = myenum.nextElement();
			out.write(contig.getDeletionReport(this.getMutantLines(), minWtCoverage, minZeroCoveragePositions));
		}
		out.close();
		
		
		
	}
	
	
	public void writeSNPReport(File outputFile, int minWtCoverage,int  minMtCoverage)throws IOException{
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));
		out.write("MutantLine\tcontig\tposition\tReferenceAllele\tAlternativeAllele\twtCoverage\tSNPCoverage\tReferenceAlleleFrequency\twtAA\tmtAA\tnonsynonymous\n");
		for(Enumeration<NLR_Contig> myenum = this.nlr_contigs.elements(); myenum.hasMoreElements();){
			NLR_Contig contig = myenum.nextElement();
			contig.filterCommonSNPs(this.getMutantLines());
			contig.filterWildtypeSNPs(this.getMutantLines(), 0.1);
			contig.filterForRegions();
			out.write(contig.getSNPReport(this.getMutantLines(), minWtCoverage, minMtCoverage));
		}
		
		out.close();
		
	}
	
	public int getNumberOfContigs(){
		return this.nlr_contigs.size();
	}
	

	public Vector<String> getMutantLines(){
		Vector<String> v = new Vector<String>();
		for(Iterator<String> iterator = mutantLines.iterator(); iterator.hasNext();){
			v.add(iterator.next());
		}
		Collections.sort(v);
		return v;
	}
	
	
	
	
	
	
	
	
	
	

	/* ******************************************************** *
	 * ********  export and import function        ************ * 
	 * ******************************************************** */
	
	
	
	/**
	 * Export the contig list to an xml file.
	 * 
	 * @param xmlFile
	 * @throws IOException
	 * @throws ParserConfigurationException
	 * @throws TransformerConfigurationException
	 * @throws TransformerException
	 */
	public void exportToXML(File xmlFile)throws IOException, ParserConfigurationException, TransformerConfigurationException, TransformerException{
		DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
		DocumentBuilder db = dbf.newDocumentBuilder();
		Document dom = db.newDocument();
		
		Element rootElement = dom.createElement("MutantHunter");
		rootElement.setAttribute("version", version+"");
		dom.appendChild(rootElement);
		
		
		Element contigsElement = dom.createElement("NBLRR_Contigs");
		rootElement.appendChild(contigsElement);
		for(Enumeration<String> myenum1 = nlr_contigs.keys(); myenum1.hasMoreElements();){
			String key = myenum1.nextElement();
			NLR_Contig contig = nlr_contigs.get(key);
			contigsElement.appendChild(contig.getXMLElement(dom));
			
			
		}
		
		
		DOMSource source = new DOMSource(rootElement) ;
		StreamResult result = new StreamResult(xmlFile);
		Transformer transformer = TransformerFactory.newInstance().newTransformer();
		if (dom.getDoctype() != null) {
		    String systemValue = (new File (dom.getDoctype().getSystemId())).getName();
		    transformer.setOutputProperty(OutputKeys.DOCTYPE_SYSTEM, systemValue);
		}
		transformer.setOutputProperty(OutputKeys.INDENT, "yes"); //without this there is no newline after elements.
		
		transformer.setOutputProperty("{http://xml.apache.org/xslt}indent-amount", "2");  //not sure what happens here but this makes the nice whitespace hierarchy in the output txt.
		transformer.transform(source, result);
		
	}
	
	
	
	/**
	 * 
	 * Read an xml file and add the data to this object. This can be used tobuild the list of NLR_Contigs or just to add data vrom additional mutant lines
	 * 
	 * @param xmlFile
	 * 			The location of the file
	 * @param addNewContigs
	 * 			if true, additional contigs that do not already exist will be added.
	 * 			if false, only new information will be added to existing NLR_Contigs.
	 * 
	 * @throws IOException
	 * @throws SAXException
	 * @throws ParserConfigurationException
	 */
	public void addXML(File xmlFile, boolean addNewContigs)throws IOException, SAXException, ParserConfigurationException{
				
		DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
		DocumentBuilder db = dbf.newDocumentBuilder();
		Document dom = db.parse(xmlFile);
		
		
		Element rootElement = dom.getDocumentElement();
		
		
		Element nlrContigsElement = (Element) rootElement.getElementsByTagName("NBLRR_Contigs").item(0);
		if( nlrContigsElement != null){
			NodeList nodeList = nlrContigsElement.getElementsByTagName("NBLRR_Contig");
			for( int i = 0; i< nodeList.getLength(); i++){
				Element nlrContigElement = (Element) nodeList.item(i);
				NLR_Contig contig = new NLR_Contig(nlrContigElement);
				if( nlr_contigs.containsKey(contig.getContigName())){
					nlr_contigs.get(contig.getContigName()).mergeContig(contig);
				}else{
					if( addNewContigs){
						nlr_contigs.put(contig.getContigName(), contig);
					}	
				}
				
				
			}
		}
		
		for(Enumeration<String> myenum = this.nlr_contigs.keys(); myenum.hasMoreElements();){
			String key = myenum.nextElement();
			for(Enumeration<String> myenum2 = this.nlr_contigs.get(key).getMutantLineNames().elements(); myenum2.hasMoreElements();){
				mutantLines.add(myenum2.nextElement());
			}
		}
		
	}
	
	
	
	
	
	
	
	
	
	

	
	
	public static void main(String[] args){
		
		CLI cli = new CLI();
		cli.parseOptions(args);
		
		
		
		
		try{
			
			if( !cli.hasOption("w") || !cli.hasOption("m") || !cli.hasOption("b") || !cli.hasOption("o")){
				throw new CLIParseException("Missing required option: -w, -m, -b and -o are required.");
			}
			
			File wtFile = new File(cli.getArg("w"));
			Vector<String> mutants = cli.getArgs("m");
			File blastFile = new File(cli.getArg("b"));
			File outputFile = new File(cli.getArg("o"));
			
			
			MutantHunter hunter = new MutantHunter(wtFile);
			for(Enumeration<String> myenum = mutants.elements(); myenum.hasMoreElements();){
				String mutant = myenum.nextElement();
				hunter.addXML(new File(mutant),false);
			}
			
			hunter.addRegions(blastFile, 100, false);
			
			
			
			int minMutants = 2;
			
			if( cli.hasOption("n")){
				try{
					minMutants = Integer.parseInt(cli.getArg("n"));
				}catch(NumberFormatException e){
					throw new CLIParseException("Argument for -n has to be an int");
				}
			}
			
			
			int minWtCov = 10;
			int minCoverageToConsiderSNP = 10;
			
			if( cli.hasOption("c")){
				try{
					minWtCov = Integer.parseInt(cli.getArg("c"));
					minCoverageToConsiderSNP = Integer.parseInt(cli.getArg("c"));
				}catch(NumberFormatException e){
					throw new CLIParseException("Argument for -c has to be an int");
				}
			}
			
					
			double maxRefAlleleFrequency = 0.1;
			if( cli.hasOption("a")){
				try{
					maxRefAlleleFrequency = Double.parseDouble(cli.getArg("a"));
				}catch(NumberFormatException e){
					throw new CLIParseException("Argument for -a has to be an float");
				}
			}
			
			int minNumberOfZeroCoveragePositions = 50;
			if( cli.hasOption("z")){
				try{
					minNumberOfZeroCoveragePositions = Integer.parseInt(cli.getArg("z"));
				}catch(NumberFormatException e){
					throw new CLIParseException("Argument for -z has to be an int");
				}
			}
			
			
			
			
			
			
			
			
			hunter.findCandidates(outputFile, minWtCov, maxRefAlleleFrequency, minCoverageToConsiderSNP, minNumberOfZeroCoveragePositions, minMutants, false);
			
			
			
			
			
			
			
			
			
			
			
		}catch (CLIParseException e){
			e.printStackTrace();
			
			String s = "-w <wt.xml>\t\t\tThe XML file generated with Pileup2XML.jar made from wildtyp\n"+
					   "-m <mt1.xml [mtn.xml]*\t\tThe XML files generated with Pileup2XML.jar made from mutants\n"+
					   "-b <blast.xml>\t\t\tThe blast file of contigs vs. baits. Use NCBI blast+ and -outfmt 5\n"+
					   "-o <output.txt>\t\t\tOutput file with andidates\n"+
					   "-n <int>\t\t\tMinimum number of mutants to report a contig. Default is 2\n"+
					   "-c <int>\t\t\tMininum coverage for mappings to be regarded. Default is 10\n"+
					   "-a <float>\t\t\tMaximum reference allele frequency to consider a SNP. Default is 0.1\n"+
					   "-z <int>\t\t\tNumber of coherent positions with zero coverage to call a deletion mutant. Default is 50\n";
			
			System.err.println(s);
			
		}
		catch (ParserConfigurationException e){
			e.printStackTrace();
		}
		catch (IOException e){
			e.printStackTrace();
		}
		catch (SAXException e){
			e.printStackTrace();
		}
		
		
		
	
	}
	
	
	
	
	
	
	
	
	
	
}
