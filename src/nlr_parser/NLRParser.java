package nlr_parser;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.concurrent.ExecutionException;

import javax.xml.parsers.ParserConfigurationException;
import javax.xml.transform.TransformerException;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import org.xml.sax.SAXException;

import supportClasses.BioSequence;
import supportClasses.FastaReader;
import supportClasses.FastqReader;

public class NLRParser {
	
	
	static final double VERSION = 2.0;
	
	Hashtable<String, MastMotifHitList> lists;
	
	
	
	
	public NLRParser(File xmlFile)throws SAXException, ParserConfigurationException, IOException{
		this.lists = MastMotifHitList.importFromXML(xmlFile);
	}
	
	public NLRParser(File inputFile, File mastExe, File memeXML, int numThreads, int numSequencesPerCall, double pvalue_threshold)throws ExecutionException, InterruptedException, IOException{
		MastParallelExecuter mpe = new MastParallelExecuter(mastExe, memeXML, numThreads, numSequencesPerCall);
		
		if(isFasta(inputFile)){
			if(isDNA(inputFile)){
				this.lists = mpe.executeMastDNA(inputFile, pvalue_threshold);
			}else{
				this.lists = mpe.executeMastProtein(inputFile, pvalue_threshold);
			}
		}else{
			if( isFastq(inputFile)){
				this.lists = mpe.executeMastDNAFastq(inputFile, pvalue_threshold);
			}
		}
		
		mpe.close();
	}
	
	
	
	
	
	
	
	
	public static double getVersion(){
		return VERSION;
	}
	
	public Hashtable<String,MastMotifHitList> getMastMotifHitLists(){
		return this.lists;
	}
	
	
	
	
	
	
	
	
	public static void main(String[] args){
		
		
		
		Options options = new Options();
		options.addOption( OptionBuilder.withLongOpt( "input")
										.withDescription("input file in fasta, fastq or xml format. If this is a sequence file, parameters y and x are required as well. If this an xml file, it is assumed that this is a backup from NLR-Parser.")
										.hasArg()
										.create('i'));
		options.addOption( OptionBuilder.withLongOpt( "output")
										.withDescription("output File")
										.hasArg()
										.withArgName("output.tsv")
										.create('o'));
		options.addOption(OptionBuilder.withLongOpt("writeGFF")
										.withDescription("Output general feature format (GFF)")
										.hasArg()
										.withArgName("output.gff")
										.create('g'));
		options.addOption(OptionBuilder.withLongOpt("backup")
									   .withDescription("A backup xml file")
									   .hasArg()
									   .withArgName("backup.xml")
									   .create('c'));
		options.addOption(OptionBuilder.withLongOpt("bed")
				   						.withDescription("Output in BED format")
				   						.hasArg()
				   						.withArgName("output.bed")
				   						.create('b'));
				options.addOption(OptionBuilder.withLongOpt("pValue")
									   .withDescription("p-value threshold for each individual motif.")
									   .hasArg()
									   .withArgName("pvalue")
									   .create('p'));
		options.addOption(OptionBuilder.withLongOpt("mast")
										.withDescription("The path to the mast program")
										.hasArg()
										.create('y'));
		options.addOption(OptionBuilder.withLongOpt("meme")
										.withDescription("The path to the meme.xml file containing the NLR motifs")
										.hasArg()
										.create('x'));
		options.addOption(OptionBuilder.withLongOpt("numThreads")
										.withDescription("Number of Threads")
										.hasArg()
										.create('t'));
		options.addOption(OptionBuilder.withLongOpt("numSequencesPerMastCall")
										.withDescription("Number of sequences that are submitted to mast analysis at the same time. Since the pvalue depends on the length of input data set you might miss hits if you sequence n is too large. Default 1000")
										.hasArg()
										.create('n'));
		options.addOption(OptionBuilder.withLongOpt("verbous")
										.withDescription("verbous. Print tsv output format to standard.out")
										.create('v'));
		options.addOption(OptionBuilder.withLongOpt("help")
				   						.withDescription("display this help")
				   						.create('h'));
		
		
		CommandLineParser parser = new PosixParser();
		//i,o,g,b,p,y,x,t,n,h
		
		try {
			
			CommandLine line = parser.parse( options, args );
			
			if( line.hasOption('h')){
				HelpFormatter formatter = new HelpFormatter();
		        formatter.printHelp( "SequenceTools", options );
			}else{
				if(!line.hasOption('i')){
					throw new ParseException ("no input file found. Provide input file with parameter -i");
				}
				File inputFile = new File(line.getOptionValue('i'));
				
				
				NLRParser nlrParser;
				
				
				boolean isXML = false;
				BufferedReader in = new BufferedReader(new FileReader(inputFile));
				String inputLine = in.readLine();
				while( inputLine!= null && inputLine.trim().length()==0){
					inputLine =in.readLine();
				}
				in.close();
				if( inputLine.toLowerCase().contains("<?xml") ){
					isXML = true;
				}
				
				
				
				
				
				if(isXML){
					nlrParser = new NLRParser(inputFile);
				}else{
					if(!line.hasOption('y') || !line.hasOption('x')){
						throw new ParseException("Input File is not XML format. If this is a sequence file you have to provide arguments -x and -y");
					}
					File mastExe = new File(line.getOptionValue('y'));
					File memeXML = new File(line.getOptionValue('x'));
					if(!memeXML.exists()){
						throw new ParseException("meme.xml not found");
					}
					if(!mastExe.exists()){
						throw new ParseException("MAST not found");
					}
					
					int numThreads = 1;
					if( line.hasOption('t')){
						try{
							numThreads = Integer.parseInt(line.getOptionValue('t'));
						}catch(NumberFormatException e){
							System.err.println("WARNING: Wrong parameter for -t. Integer was expected. Program executed with -t 1.");
						}
					}
					
					double pvalue = 1E-5;
					if( line.hasOption('p')){
						try {
							pvalue = Double.parseDouble(line.getOptionValue('p'));
						} catch (Exception e) {
							System.err.println("WARNING: Wrong parameter for -p. Float was expected. Program executed with -p 1e-5.");
						}
						
					}
					
					
					
					int numberOfSeqeuncesPerMastCall = 1000;
					if(line.hasOption('n')){
						try{
							numberOfSeqeuncesPerMastCall = Integer.parseInt(line.getOptionValue('n'));
						}catch ( NumberFormatException e){
							System.err.println("WARNING: Wrong parameter for -n. Integer was expected. Program executed with -n 1000.");
						}
					}
					
					
					nlrParser = new NLRParser(inputFile, mastExe, memeXML, numThreads, numberOfSeqeuncesPerMastCall, pvalue);
					
					nlrParser.addSequence(inputFile);
					
					
				}
				
				
				
				
				if( line.hasOption('o')){
					
					File outputFile = new File(line.getOptionValue('o'));
					if( outputFile.exists()){
						outputFile.renameTo(new File(outputFile.getParentFile(), outputFile.getName()+".bak"));
						outputFile = new File(line.getOptionValue('o'));
						
						
						File intermediateFile = new File(outputFile.getAbsolutePath());
						while(intermediateFile.exists()){
							intermediateFile = new File(intermediateFile.getAbsolutePath() + ".bak");
						}
						if(outputFile.exists()){
							outputFile.renameTo(intermediateFile);
							System.err.println("WARNING: Output file exists. The existing file was renamed to "+intermediateFile.getName());
						}
						
						outputFile = new File(line.getOptionValue('o'));
					}
					nlrParser.writeTSV(outputFile);
				}
				
				if(line.hasOption('g')){
					
					File outputFile = new File(line.getOptionValue('g'));
					if( outputFile.exists()){
						outputFile.renameTo(new File(outputFile.getParentFile(), outputFile.getName()+".bak"));
						outputFile = new File(line.getOptionValue('g'));
						
						
						File intermediateFile = new File(outputFile.getAbsolutePath());
						while(intermediateFile.exists()){
							intermediateFile = new File(intermediateFile.getAbsolutePath() + ".bak");
						}
						if(outputFile.exists()){
							outputFile.renameTo(intermediateFile);
							System.err.println("WARNING: Output file exists. The existing file was renamed to "+intermediateFile.getName());
						}
						
						outputFile = new File(line.getOptionValue('g'));
					
				}
					nlrParser.writeGFF(outputFile);
				}
				
				if( line.hasOption('v')){
					nlrParser.writeStdOut();
				}
				
				if(line.hasOption('c')){
					nlrParser.writeBackup(new File(line.getOptionValue('c')));
				}
				
				
				
				
				
				if(line.hasOption('b')){
					
					File outputFile = new File(line.getOptionValue('b'));
					if( outputFile.exists()){
						outputFile.renameTo(new File(outputFile.getParentFile(), outputFile.getName()+".bak"));
						outputFile = new File(line.getOptionValue('b'));
						
						
						File intermediateFile = new File(outputFile.getAbsolutePath());
						while(intermediateFile.exists()){
							intermediateFile = new File(intermediateFile.getAbsolutePath() + ".bak");
						}
						if(outputFile.exists()){
							outputFile.renameTo(intermediateFile);
							System.err.println("WARNING: Output file exists. The existing file was renamed to "+intermediateFile.getName());
						}
						
						outputFile = new File(line.getOptionValue('b'));
					
				}
					nlrParser.writeBED(outputFile);
				}
				
				
				
				
				
			}
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
					
			
			
			
			
			
			
			
		} catch (ParseException   e) {
			e.printStackTrace();
			
			HelpFormatter formatter = new HelpFormatter();
	        formatter.printHelp( "SequenceTools", options );
			
		}
		
		catch (IOException   e) {
			e.printStackTrace();
		}
		
		catch (ExecutionException   e) {
			e.printStackTrace();
		}
		
		catch (InterruptedException   e) {
			e.printStackTrace();
		}
		catch(ParserConfigurationException e){
			e.printStackTrace();
		}
		catch (SAXException e){
			e.printStackTrace();
		}
		catch(TransformerException e){
			
		}
	}
	
	
	
	
	public static Hashtable<String,MastMotifHitList> executeNLRParser(File inputFile, File mastExe, File memeXml, int numThreads, int numberOfSeqeuncesPerMastCall, boolean isFastq, boolean isProtein, double pvalue  ) throws ExecutionException,InterruptedException,IOException{
		
		
		
		MastParallelExecuter mpe = new MastParallelExecuter( mastExe,  memeXml,  numThreads,  numberOfSeqeuncesPerMastCall);
		System.err.println("DEBUG: Temporary directory is " + mpe.getTempDir().getAbsolutePath());
		
		Hashtable<String,MastMotifHitList> resultSet;
		
		if( isFastq ){
			System.err.println("DEBUG: dna fastq");
			resultSet =  mpe.executeMastDNAFastq(inputFile, pvalue);
		}else{
			if( isProtein){
				System.err.println("DEBUG: protein fasta");
				resultSet = mpe.executeMastProtein(inputFile, pvalue);
			}else{
				System.err.println("DEBUG: dna fasta");
				resultSet = mpe.executeMastDNA(inputFile, pvalue);
			}
			
			
		}
		System.err.println("DEBUG: Size of the result set: " + resultSet.size());
		
		mpe.close();
		return resultSet;
	}
	
	
	
	public  Hashtable<String,MastMotifHitList> addSequence( File sequenceFile)throws IOException{
			
		
		if(isFasta(sequenceFile)){
			
			FastaReader fastaReader = new FastaReader(sequenceFile);
			for (BioSequence seq = fastaReader.readEntry(); seq != null; seq = fastaReader.readEntry()) {
				if(seq.isDNA()){
					if( this.lists.containsKey(seq.getIdentifier() + "+")){
						this.lists.get(seq.getIdentifier() + "+").addDnaSequence(seq.getSequence());
						
					}
					if( this.lists.containsKey(seq.getIdentifier() + "-")){
						this.lists.get(seq.getIdentifier()+ "-").addDnaSequence(seq.getReverseComplementarySequence());
					}
					
				}else{
					this.lists.get(seq.getIdentifier()).addAASequence(seq.getSequence());
						
				}
			}
			fastaReader.close();
		}else if(isFastq(sequenceFile)){
			FastqReader fastaReader = new FastqReader(sequenceFile);
			for (BioSequence seq = fastaReader.readEntry(); seq != null; seq = fastaReader.readEntry()) {
				if( this.lists.containsKey(seq.getIdentifier() + "+")){
					this.lists.get(seq.getIdentifier() + "+").addDnaSequence(seq.getSequence());
				}
				if( this.lists.containsKey(seq.getIdentifier() + "-")){
					this.lists.get(seq.getIdentifier()+ "-").addDnaSequence(seq.getReverseComplementarySequence());
				}
			}
			fastaReader.close();
		}
		
		return this.lists;
	}
	
	
	private static boolean isDNA(File inputFile)throws IOException{
		boolean isDNA = false;
		FastaReader reader = new FastaReader(inputFile);
		BioSequence seq = reader.readEntry();
		if(seq.isDNA()){
			isDNA = true;
		}
		reader.close();
		return isDNA;
	}
	
	
	
	
	private static boolean isFasta(File inputFile)throws IOException {
		BufferedReader in = new BufferedReader(new FileReader(inputFile));
		
		for (String inputline = in.readLine();inputline != null;inputline = in.readLine() ) {
			if(inputline.startsWith("#")){
				continue;
			}
			if(inputline.trim().equalsIgnoreCase("")){
				continue;
			}
			if(inputline.trim().startsWith(">")){
				in.close();
				return true;
			}else{
				in.close();
				return false;
			}
		}
		in.close();
		
		return false;
	}
	
	
	/**
	 * check if the first line starts with an @ and the third with a + if that's the case we assume Fastq. 
	 * @param inputFile
	 * @return
	 * @throws IOException
	 */
	private static boolean isFastq(File inputFile)throws IOException {
		BufferedReader in = new BufferedReader(new FileReader(inputFile));
		
		for (String inputline = in.readLine();inputline != null;inputline = in.readLine() ) {
			
			if(inputline.trim().equalsIgnoreCase("")){
				continue;
			}
			if(inputline.trim().startsWith("@")){
				in.readLine();
				inputline = in.readLine();
				if(inputline.startsWith("+")){
					in.close();
					return true;
				}else{
					in.close();
					return false;
				}
				
			}else{
				in.close();
				return false;
			}
		}
		in.close();
		
		return false;
	}
	
	public  void writeTSV( File outputFile)throws IOException{
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));
		out.write("SequenceName\tclass\tcomplete\tstrand\tstart\tend\t1st-motif-frame\tNB-ARC\t2-NBLRR-Signal\tMotifList");
		out.newLine();
		for(Enumeration<String> myenum = this.lists.keys(); myenum.hasMoreElements();){
			String key = myenum.nextElement();
			MastMotifHitList list = this.lists.get(key);
			
			if( list.isNibbler()){
				
				String sequenceName = key.substring(0, key.length()-1);	
				String strand = "N/A";
				String start = list.getAAStart() + "";
				String end   = list.getAAEnd() + "";
				
				
				if(key.substring( key.length()-1).equalsIgnoreCase("-")){
					strand = "reverse";
					start = "N/A";
					end   = "N/A";
					if( list.getDNASequence() != null && list.getDNASequence().length() > 0){
						int length = list.getDNASequence().length();
						start = ( length - list.getDNAEnd()) + "";
						end   = (length - list.getDNAStart()) +"";
					}
					
					
				}else if(key.substring( key.length()-1).equalsIgnoreCase("+")){
					strand = "forward";
					start = list.getDNAStart() + "";
					end   = list.getDNAEnd() + "";
				}
				
				
				
				
				
				
				String nlrclass = list.getNBLRR_Class();
				String nbarc = list.getNBARC_Sequence();
				
				String isComplete = "partial";
				if(list.isComplete()){
					isComplete = "complete";
				}
				
				out.write( sequenceName + "\t" + nlrclass +"\t" + isComplete + "\t"+ strand+"\t" + start + "\t" + end + "\t" + list.getFirstMotifFrame() + "\t" + nbarc +"\t"+list.hasPotentiallyTwoNBLRRs() +"\t"+list.getMotifListString());
				out.newLine();
				
				
			}
			
		}
		
		
		out.close();
	}
	
	
	
	private  void writeStdOut( )throws IOException{
		
		System.out.println("SequenceName\tclass\tcomplete\tstrand\tstart\tend\t1st-motif-frame\tNB-ARC\t2-NBLRR-Signal\tMotifList");
		
		for(Enumeration<String> myenum = this.lists.keys(); myenum.hasMoreElements();){
			String key = myenum.nextElement();
			MastMotifHitList list = this.lists.get(key);
			//System.out.println(key);
			if( list.isNibbler()){
				
				String sequenceName = key.substring(0, key.length()-1);	
			
				String strand = "N/A";
				String start = list.getAAStart() + "";
				String end   = list.getAAEnd() + "";
				
				
				if(key.substring( key.length()-1).equalsIgnoreCase("-")){
					strand = "reverse";
					start = "N/A";
					end   = "N/A";
					if( list.getDNASequence() != null && list.getDNASequence().length() > 0){
						int length = list.getDNASequence().length();
						start = ( length - list.getDNAEnd()) + "";
						end   = (length - list.getDNAStart()) +"";
					}
					
					
				}else if(key.substring( key.length()-1).equalsIgnoreCase("+")){
					strand = "forward";
					start = list.getDNAStart() + "";
					end   = list.getDNAEnd() + "";
				}
				
				
				
				
				
				String nlrclass = list.getNBLRR_Class();
				String nbarc = list.getNBARC_Sequence();
				
				String isComplete = "partial";
				if(list.isComplete()){
					isComplete = "complete";
				}
				
				System.out.println( sequenceName + "\t" + nlrclass +"\t" + isComplete + "\t"+ strand+"\t" + start + "\t" + end + "\t" + list.getFirstMotifFrame() + "\t" + nbarc +"\t"+list.hasPotentiallyTwoNBLRRs() +"\t"+list.getMotifListString());
				
				
				
			}
			
		}
		
		
		
	}
	
	
	
	
	/**
	 * This is experimental. Output GFF.
	 * 
	 * 
	 * 
	 * @param resultSet
	 * @param outputFile
	 * @throws IOException
	 */
	public  void writeGFF( File outputFile)throws IOException{
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));

		out.write("##gff-version 2");
		out.newLine();
		out.write("##source-version NLR-Parser " + NLRParser.getVersion());
		out.newLine();
		
		DateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd HH:mm");
		Date date = new Date();
		out.write("##date "+dateFormat.format(date)); 
		out.newLine();
		out.write("##Type DNA");
		out.newLine();
		out.write("#seqname\tsource\tfeature\tstart\tend\tscore\tstrand\tframe\tattribute");
		out.newLine();
		
		
		for(Enumeration<String> myenum1 = this.lists.keys(); myenum1.hasMoreElements();){
			String key = myenum1.nextElement();
			String sequenceName = key.substring(0, key.length()-1);
			MastMotifHitList list = this.lists.get(key);
			int length= 0;
			if( list.getDNASequence()!= null){
				length = list.getDNASequence().length();
			}
			
			for( Enumeration<MastMotifHit> myenum2 = list.getMotifs().elements(); myenum2.hasMoreElements();){
				MastMotifHit hit = myenum2.nextElement();
				out.write(sequenceName +"\tNLR-Parser\tMastMotif\t"  );
				
				int frame = hit.getFrame();
				
				
				boolean isForwardStrand = hit.forwardStrand;
				
				int nucl_start = hit.getDNAPos();
				int nucl_end = hit.getDNAPos() + hit.getMatch().length()*3;
				String strand = "+";
				
				if(!isForwardStrand){
					strand = "-";
					nucl_start = length - nucl_end;
					nucl_end = length - hit.getDNAPos();
				}
				
				
				out.write(nucl_start + "\t" + nucl_end + "\t" +  hit.getPvalue() + "\t" +strand +"\t" + (Math.abs(frame) ) + "\t" + "name "+hit.getMotif() );
				out.newLine();
				
				
				
			}
			
		}	
		
		
		
		out.close();
		
		
		
	}
	
	/*
	 public void exportToBED(File outputFile)throws IOException{
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));
		out.write("track name=\"NLR_Motifs\"");
		out.newLine();
		out.write("itemRgb=\"On\"");
		out.newLine();
		
		for(Enumeration<String> myenum = this.motifLists.keys(); myenum.hasMoreElements();){
			String key = myenum.nextElement();
			NLR_MotifList list = this.motifLists.get(key);
			out.write(list.getBED());
		}
		
		out.close();
	 */
	
	public void writeBED(File outputFile)throws IOException{
		
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));
		out.write("#track name=\"NLR_Motifs\"");
		out.newLine();
		out.write("#itemRgb=\"On\"");
		out.newLine();
		
		for(Enumeration<String> myenum1 = this.lists.keys(); myenum1.hasMoreElements();){
			String key = myenum1.nextElement();
			String sequenceName = key.substring(0, key.length()-1);
			MastMotifHitList list = this.lists.get(key);
			int length= 0;
			if( list.getDNASequence()!= null){
				length = list.getDNASequence().length();
			}
			
			for( Enumeration<MastMotifHit> myenum2 = list.getMotifs().elements(); myenum2.hasMoreElements();){
				MastMotifHit hit = myenum2.nextElement();
				
				
				boolean isForwardStrand = hit.forwardStrand;
				
				int nucl_start = hit.getDNAPos();
				int nucl_end = hit.getDNAPos() + hit.getMatch().length()*3;
				String strand = "+";
				
				if(!isForwardStrand){
					strand = "-";
					nucl_start = length - nucl_end;
					nucl_end = length - hit.getDNAPos();
				}
				
				int[] a = hit.getRGBMotifColor();
				out.write(sequenceName + "\t" + nucl_start + "\t" + nucl_end + "\t" + hit.getMotif() + "\t0\t" + strand + "\t" + nucl_start + "\t" + nucl_end +"\t" + a[0]+","+a[1]+","+a[2] );
				
				out.newLine();
				
			}
			
		}
		
		
		out.close();
		
	}
	
	
	public void writeBackup(File outputFile)throws TransformerException, ParserConfigurationException, IOException{
		MastMotifHitList.exportToXML(this.lists, outputFile);
	}
	
	
	
	/**
	 * Deliver the color code for motifs as we know it from the mast html output. 
	 * 
	 * @param motifNumber
	 * @return
	 */
	public static int[] getRGBMotifColor(int motifNumber){
		if( motifNumber ==1 ){ int[] a = {0,255,255}; return a;}
		if( motifNumber ==2 ){ int[] a = {0,0,255}; return a;}
		if( motifNumber ==3 ){ int[] a = {255,0,0}; return a;}
		if( motifNumber ==4 ){ int[] a = {255,0,255}; return a;}
		if( motifNumber ==5 ){ int[] a = {255,255,0}; return a;}
		if( motifNumber ==6 ){ int[] a = {0,255,0}; return a;}
		if( motifNumber ==7 ){ int[] a = {0,128,128}; return a;}
		if( motifNumber ==8 ){ int[] a = {68,68,68}; return a;}
		if( motifNumber ==9 ){ int[] a = {0,128,0}; return a;}
		if( motifNumber ==10 ){ int[] a = {192,192,192}; return a;}
		if( motifNumber ==11 ){ int[] a = {128,0,128}; return a;}
		if( motifNumber ==12 ){ int[] a = {128,128,0}; return a;}
		if( motifNumber ==13 ){ int[] a = {0,0,128}; return a;}
		if( motifNumber ==14 ){ int[] a = {128,0,0}; return a;}
		if( motifNumber ==15 ){ int[] a = {255,255,255}; return a;}
		if( motifNumber ==16 ){ int[] a = {0,255,255}; return a;}
		if( motifNumber ==17 ){ int[] a = {0,0,255}; return a;}
		if( motifNumber ==18 ){ int[] a = {255,0,0}; return a;}
		if( motifNumber ==19 ){ int[] a = {255,0,255}; return a;}
		if( motifNumber ==20 ){ int[] a = {255,255,0}; return a;}
		
		return new int[3];
	}
	

}
