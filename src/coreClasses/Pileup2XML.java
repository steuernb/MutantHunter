package coreClasses;

import java.io.File;
import java.io.IOException;

import javax.xml.parsers.ParserConfigurationException;
import javax.xml.transform.TransformerException;

import nlr_parser.NLRParser;

import org.xml.sax.SAXException;

import supportClasses.CLI;
import supportClasses.CLIParseException;

public class Pileup2XML {

	

	
	
	public static void main(String[] args){
		
		
		CLI cli = new CLI();
		cli.parseOptions(args);
		
		String helpString =		"-i <inputPileup>\t\t\tInput file in mpileup format\n"+
								"-m <nlrparser.xml>\t\t\tThe annotation file of the assembly. Run NLR-Parser with option -c\n"+
								"-o <output.xml>\t\t\tThe output file in xml format\n"+
								"-w\t\t\tThis pileup belongs to the wildtype\n"+
								"-c <int>\t\t\tminimum coverage to report a SNP\n";
								

		
		try{
			
			if( !cli.hasOption("i") || !cli.hasOption("m")){
				throw new CLIParseException("Missing option. -i and -m are mandatory.");
			}
			
			
			File inputPileupFile = new File(cli.getArg("i"));
			
			File outputXMLFile = new File( inputPileupFile.getParentFile(), inputPileupFile.getName().split("\\.pileup")[0]+".xml");
			
			if(cli.hasOption("o")){
				outputXMLFile = new File(cli.getArg("o"));
			}
			
			
			boolean isWildtype = false;
			if(cli.hasOption("w")){
				isWildtype = true;
			}
			
			int minCoverage = 5;
			if(cli.hasOption("c")){
				try{
					minCoverage = Integer.parseInt(cli.getArg("c"));
				}catch (NumberFormatException e){System.err.println("Warning: wrong input value for -c. This needs to be in int. Default is used.");};	
			}
			
			
			MutantHunter hunter;
			//option1 There is a mast backup file
			
			NLRParser nlrparser = new NLRParser(new File(cli.getArg("m")));
				
			hunter = new MutantHunter();
			hunter.setNLRContigList(nlrparser.getMastMotifHitLists());
						
			
			
			String mutantLine = inputPileupFile.getName();
			if( isWildtype){
				mutantLine = MutantHunter.wildtype;
			}
			hunter.readPileupFile(mutantLine, inputPileupFile, minCoverage);
			
			
			hunter.exportToXML(outputXMLFile);
			
		}catch(IOException e){
			e.printStackTrace();
			System.err.println(helpString);
		}
		
		catch(SAXException e){
			e.printStackTrace();
		}
		catch(ParserConfigurationException e){
			e.printStackTrace();
		}
		
		catch (TransformerException e){
			e.printStackTrace();
		}
		catch(CLIParseException e){
			e.printStackTrace();
			System.err.println(helpString);
		}
	}
	
	
	
	
	
	
	
	
	
	
	
	
}
