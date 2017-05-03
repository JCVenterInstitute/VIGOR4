package com.vigor.test;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import org.apache.commons.configuration2.INIConfiguration;
import org.apache.commons.configuration2.builder.fluent.Configurations;
import org.springframework.core.io.ClassPathResource;
import org.springframework.core.io.Resource;

import com.vigor.utils.VigorUtils;

public class AppTest {
	
	public static void main(String[] args) {
		Map<String,String> vigorParametersList = new HashMap<String,String>();
		Map<String,String> map = new HashMap<String,String>();
		try{
			
			Configurations configs = new Configurations();
			
			Resource resource = new ClassPathResource(VigorUtils.getConfigIniPath());
			
			INIConfiguration iniConfig = configs.ini(resource.getFile());
			
			iniConfig.getSections();
			
			Set<String> setOfSections = iniConfig.getSections();

			/*Iterator sectionNames = setOfSections.iterator();

			while (sectionNames.hasNext()) {
				String sectionName = sectionNames.next().toString();
				SubnodeConfiguration sObj = iniConfig.getSection(sectionName);
				Iterator it1 = sObj.getKeys();
				
				while (it1.hasNext()) {
					// Get element
					Object key = it1.next();
					vigorParametersList.put(key.toString(), sObj.getString(key.toString()));
				}

			}
			
			for (String key : vigorParametersList.keySet()) {
				System.out.println(key+" : "+vigorParametersList.get(key));
			}*/
			
			
			setOfSections.stream().forEach(i -> iniConfig.getSection(i).getKeys().forEachRemaining(n -> map.put(n, iniConfig.getSection(i).getString(n))));
			System.out.println("*********************************");
			for (String key : map.keySet()) {
				System.out.println(key+" : "+vigorParametersList.get(key));
			}
			
			
			
		}catch (Exception e) {
			e.printStackTrace();
		}
	}

}
