package org.jcvi.vigor.utils;

public class OSValidator {
	
	private static String OS = System.getProperty("os.name").toLowerCase();

	public static Boolean is64Bit() {
		String arch = System.getenv("PROCESSOR_ARCHITECTURE");
		String wow64Arch = System.getenv("PROCESSOR_ARCHITEW6432");

		return arch.endsWith("64") || wow64Arch != null && wow64Arch.endsWith("64") ? true : false;
	}

	public static boolean isWindows() {
		return (OS.indexOf("win") >= 0);
	}

	public static boolean isMac() {
		return (OS.indexOf("mac") >= 0);
	}

	public static boolean isUnix() {
		return (OS.indexOf("nix") >= 0 || OS.indexOf("nux") >= 0 || OS.indexOf("aix") > 0);
	}

	public static boolean isSolaris() {
		return (OS.indexOf("sunos") >= 0);
	}

}
