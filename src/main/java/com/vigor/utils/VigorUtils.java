package com.vigor.utils;

import java.io.File;
import java.util.regex.Pattern;

public class VigorUtils {

	private static String BLAST_CLASSPATH = "vigorResources" + File.separator + "blast" + File.separator;
	private static String WINDOWS = "windows";
	private static String LINUX = "linux";

	public static boolean isNullOrEmpty(Object value) {
		if (value != null && !value.toString().trim().equalsIgnoreCase(""))
			return false;
		return true;
	}

	public static String removeSpaces(String value) {
		if (value != null && !value.toString().trim().equalsIgnoreCase("")) {
			return value.replaceAll("(?m)^[ \t]*\r?\n", "").trim();
		}
		return "";
	}

	public static Boolean isDefLine(String value) {
		if (value != null && !value.toString().trim().equalsIgnoreCase("")) {
			return Pattern.compile("^>.*").matcher(value.trim()).matches();
		}
		return false;
	}

	public static Boolean isNucleotide(String value) {
		if (value != null && !value.toString().trim().equalsIgnoreCase("")) {
			return Pattern.compile("^[a-zA-Z]*$").matcher(value.trim()).matches();
		}
		return false;
	}

	
	public static String getVigorWorkSpace() {
		File theDir = new File("VigorWorkSpace");
		if (!theDir.exists()) {
			try {
				theDir.mkdir();
			} catch (SecurityException se) {

			}
		}

		return theDir.getPath();
	}

	public static String getBlastCommand(String blastFilePath, String inputFilePath, String reference_db,
			String outputFilePath, String logFilePath) {
		return blastFilePath + " -p blastx -i " + inputFilePath + " -d " + reference_db
				+ " -e 1e-5 -M BLOSUM45 -g F -F \"\" -z 3000000 -v 0 -b 100 -m 7";
	}

	public static String getVirusDatabasePath() {
		return "vigorResources" + File.separator + "data3";
	}

	public static String getConfigIniPath() {

		String vigorIniPath = "vigorResources" + File.separator + "config" + File.separator + "vigor.ini";
		return vigorIniPath;
	}

	public static String getBlastFilePath() throws IllegalArgumentException {
		if (OSValidator.isWindows()) {
			if (OSValidator.is64Bit()) {
				return BLAST_CLASSPATH + WINDOWS + File.separator + "64" + File.separator + "blastall";
			} else {
				return BLAST_CLASSPATH + WINDOWS + File.separator + "32" + File.separator + "blastall";
			}
		} else if (OSValidator.isUnix()) {
			if (OSValidator.is64Bit()) {
				return BLAST_CLASSPATH + LINUX + File.separator + "64" + File.separator + "blastall";
			} else {
				return BLAST_CLASSPATH + LINUX + File.separator + "32" + File.separator + "blastall";
			}
		} else if (OSValidator.isMac()) {
			/*
			 * if (OSValidator.is64Bit()) { return
			 * "blast"+File.pathSeparator+"64"+File.separator+"blastall"; } else
			 * { return
			 * "blast"+File.pathSeparator+"64"+File.separator+"blastall"; }
			 */
			throw new IllegalArgumentException("No Blast File Found for Mac OS");
		} else if (OSValidator.isSolaris()) {
			/*
			 * if (OSValidator.is64Bit()) { return
			 * "blast"+File.pathSeparator+"64"+File.pathSeparator+"blastall"; }
			 * else { return
			 * "blast"+File.pathSeparator+"64"+File.pathSeparator+"blastall"; }
			 */
			throw new IllegalArgumentException("No Blast File Found for Solaris OS");
		} else {
			System.out.println("Your OS is not support!!");
			throw new IllegalArgumentException("Your OS is not support!!");
		}
	}

}
