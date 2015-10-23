package supportClasses;

import java.util.Hashtable;
import java.util.Vector;

public class CLI {

	Hashtable<String,Vector<String>> options;
	
	
	
	public CLI(){
		options = new Hashtable<String,Vector<String>>();
	}
	
	
	public boolean hasOption(String option){
		if(options.containsKey(option)){
			return true;
		}
		return false;
	}
	
	
	
	public String getArg(String option){
		return options.get(option).get(0);
	}
	
	public Vector<String> getArgs(String option){
		return options.get(option);
	}
	
	public void parseOptions(String[] args){
		
		
		boolean inOption = false;
		String optionValue = "";
		Vector<String> arguments = new Vector<String>();
		
		for( int i = 0; i< args.length; i++){
			if(args[i].startsWith("-")){
				if(inOption){
					options.put(optionValue, arguments);
				}
				inOption = false;
				optionValue = args[i];
				while(optionValue.startsWith("-")){
					optionValue = optionValue.substring(1);
				}
				arguments = new Vector<String>();
			}else if(!args[i].trim().equalsIgnoreCase("")){
				arguments.add(args[i]);
			}
		}
	}
	
	
	
	
	
}
