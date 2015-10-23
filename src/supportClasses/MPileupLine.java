package supportClasses;

import java.util.Collections;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class MPileupLine {
	String line;
	
	Hashtable<String,Integer> alleles;
	
	public MPileupLine (String line){
		this.line = line;
	}
	
	
	public String getChromosome(){
		return line.split("\t")[0];
	}
	public char getRefBase(){
		return line.split("\t")[2].toCharArray()[0];
	}
	public int getPosition(){
		return Integer.parseInt(line.split("\t")[1]);
	}
	public int getCoverage(){
		return Integer.parseInt(line.split("\t")[3]);
	}
	
	
	public Hashtable<String,Integer> getAlternativeAlleles(){
		Hashtable<String,Integer> h =new Hashtable<String,Integer> ();
		Hashtable<String,Integer> alleles = this.getAlleles();
		for(Enumeration<String> myenum = alleles.keys(); myenum.hasMoreElements();){
			String key = myenum.nextElement();
			if(!key.equalsIgnoreCase(""+this.getRefBase())){
				h.put(key,alleles.get(key) );
			}
		}
		return h;
	}
	
	/**
	 * returns the original String representation of the line in the MPileup file, without the line terminator.
	 * 
	 */
	public String toString(){
		return this.line;
	}
	
	
		
	public Hashtable<String,Integer>  getAlleles()throws IllegalStateException{
		if(this.alleles != null){
			return alleles;
		}
		
		this.alleles = new Hashtable<String,Integer>();
		String s = "";
		try{
			 s = line.split("\t")[4];
		}catch(ArrayIndexOutOfBoundsException e){return alleles;}
		char[] c = s.toCharArray();
		
		
		
		int position = 0;
		while(position<c.length){
			try {
				
			
				String allele = "";
				
				//System.out.print(position + "\t" + c[position] + "\t");
				
				if      ( c[position]== '.' || c[position] == ',' ){ allele = this.getRefBase()+"";}
				else if ( c[position]== 'A' || c[position] == 'a' ){ allele = 'A' + "";}
				else if ( c[position]== 'T' || c[position] == 't' ){ allele = 'T' + "";}
				else if ( c[position]== 'G' || c[position] == 'g' ){ allele = 'G' + "";}
				else if ( c[position]== 'C' || c[position] == 'c' ){ allele = 'C' + "";}
				else if ( c[position]== '*' ){ allele = "*" ;}      //this one is new it is a deletion in the read
				else if ( c[position]== '$'                       ){ }  //end of read symbol; do nothing
				else if ( c[position]== '^'                       ){ position++;}// beginning of read symbol. Skip the next symbol this is read mapping quality
				else if ( c[position]== '+'){
					position ++; 
					//Pattern p = Pattern.compile("([\\d]+)([AaTtGgCcNnMmYyKkWwRrSsDd]+)");
					Pattern p = Pattern.compile("([\\d]+)([\\w]+)");
					Matcher m = p.matcher(s.substring(position));
					m.find();
					String mydigit = m.group(1);
					int digit = Integer.parseInt(mydigit);
					allele = "+" + mydigit + s.substring(position+mydigit.length(), position + mydigit.length() + digit);
					position = position + mydigit.length() + digit -1;
					//System.out.println(allele);
				}else if( c[position]== '-'){
					position ++; 
					Pattern p = Pattern.compile("([\\d]+)([\\w]+)");
					Matcher m = p.matcher(s.substring(position));
					m.find();
					String mydigit = m.group(1);
					int digit = Integer.parseInt(mydigit);
					allele = "-" + mydigit + s.substring(position+mydigit.length(), position + mydigit.length() + digit);
					position = position + mydigit.length() + digit -1;
				}
					
				int num = 0;
				if( alleles.containsKey(allele)){
					num = alleles.get(allele);
				}
				num++;
				if(!allele.equalsIgnoreCase("")){
					alleles.put(allele, num);
				}
				
				
			//	System.out.println(allele +"\t" + position);
				
				position++;
			} catch (IllegalStateException e) {
				System.out.println("\n\n\n\n"+line);
				throw new IllegalStateException (e.getMessage());
			}
		}
		
		
		return this.alleles;
	}
	
	
	/**
	 * @deprecated
	 * @return
	 */
	public Hashtable<String,Integer> getAllelesOLD(){
		Hashtable<String,Integer> h = new Hashtable<String, Integer> ();
		h.put("A", 0);
		h.put("T", 0);
		h.put("G", 0);
		h.put("C", 0);
		char[] c = line.split("\t")[4].toCharArray();
		int skipNext=0;
		for( int i = 0; i< c.length; i++){
			if(skipNext>0){
				skipNext--;
			}else{
				if( c[i] == '^'){
					skipNext++;
				}
				if( c[i] == '-'){
					int skipLength = 0;
					boolean isNumber = true;
					String skipLengthString = "";
					for( int myindex = 1; isNumber; myindex++){
						skipLengthString = skipLengthString + c[i+myindex];
						try {
							skipLength = Integer.parseInt(skipLengthString);
						} catch (NumberFormatException e) {
							isNumber = false;
						}
					}
					int charsForNumber = new String(skipLength + "").length();
					skipNext = charsForNumber + skipLength;
					
				}
				
				if(c[i] =='.' || c[i] == ','){
					int a = h.get(this.getRefBase()+"");
					a++;
					h.put(this.getRefBase() + "", a);
				}
				
				
				if(c[i] =='A' || c[i] == 'a'){
					int a = h.get("A");
					a++;
					h.put("A", a);
				}
				if(c[i] =='T' || c[i] == 't'){
					int a = h.get("T");
					a++;
					h.put("T", a);
				}
				if(c[i] =='G' || c[i] == 'g'){
					int a = h.get("G");
					a++;
					h.put("G", a);
				}
				if(c[i] =='C' || c[i] == 'c'){
					int a = h.get("C");
					a++;
					h.put("C", a);
				}
				
				//implement other skipping cases
			}
			
			
		}
		
		Vector<String> remove = new Vector<String> ();
		for(Enumeration<String> myenum = h.keys(); myenum.hasMoreElements();){
			String key = myenum.nextElement();
			if(h.get(key) == 0){
				remove.add(key);
			}
		}
		for(Enumeration<String>myenum = remove.elements(); myenum.hasMoreElements();){
			h.remove(myenum.nextElement());
		}
		
		return h;
	}
	
	
	
	public boolean hasInDels(int minNumber){
		for(Enumeration<String> myenum = this.getAlleles().keys(); myenum.hasMoreElements();){
			String allele = myenum.nextElement();
			if( allele.startsWith("+") || allele.startsWith("-")){
				if(alleles.get(allele).intValue()>=minNumber){
					return true;
				}
			}
		}
		return false;
	}
	
	public boolean hasInDels(double alleleFrequency){
		for(Enumeration<String> myenum = this.getAlleles().keys(); myenum.hasMoreElements();){
			String allele = myenum.nextElement();
			if( allele.startsWith("+") || allele.startsWith("-")){
				if(alleles.get(allele).doubleValue()/(double)this.getCoverage()>=alleleFrequency){
					return true;
				}
			}
		}
		return false;
	}
	
	public boolean hasInDels(){
		return hasInDels(0);
	}
	
	public String getLargestAlternativeAllele(){
		Hashtable<String,Integer>  h =  this.getAlternativeAlleles();
		String allele = "";
		int i = 0;
		for(Enumeration<String> myenum = h.keys(); myenum.hasMoreElements();){
			String key = myenum.nextElement();
			int j = h.get(key);
			if(j>i){
				i = j;
				allele = key;
			}
		}
		
		return allele;
		
	}
	
	public int getLargestAlternativeAlleleCount(){
		Hashtable<String,Integer>  h =  this.getAlternativeAlleles();
		//String allele = "";
		int i = 0;
		for(Enumeration<String> myenum = h.keys(); myenum.hasMoreElements();){
			String key = myenum.nextElement();
			int j = h.get(key);
			if(j>i){
				i = j;
				//allele = key;
			}
		}
		
		return i;
	}
	
	
	public double getLargestAlternativeAlleleFrequency(){
		double cov = this.getCoverage();
		try{
			double allele = this.getAlternativeAlleles().get(getLargestAlternativeAllele());
			return allele / cov;
		}catch(NullPointerException e){return 0.0;}
	}
	
	public double getReferenceAlleleFrequency(){
		int ref = 0;
		if(this.getAlleles().containsKey(this.getRefBase()+"")){
			ref = this.getAlleles().get(this.getRefBase()+"");
		}
		
		return ((double) ref) / ((double)this.getCoverage()) ;
	}
	
	public Vector<Integer> getAlternativeAlleleFrequencies(){
		Vector<Integer>  v = new Vector<Integer>();
		for(Enumeration<Integer> myenum = this.getAlternativeAlleles().elements(); myenum.hasMoreElements();){
			v.add(myenum.nextElement());
		}
		Collections.sort(v);
		return v;
	}
	
}
