// Simple java program to write zlib format data as found in 
// Amiramesh files using the HxZip encoding
// (c) Copyright 2010 Gregory Jefferis. All Rights Reserved.

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.util.zip.DeflaterOutputStream;

public class WriteHxZipdata {
	
	static boolean verbose = false;
	
	public static void main(String [ ] args)
	{
		if(args.length!=2) {
			System.err.println("usage: java WriteHxZipdata <infile> <outfile>");
			return;
		}

		String infilepath=args[0];
		String outfilepath=args[1];
		File infile = new File(infilepath);
		File outfile = new File(outfilepath);

		try{
			
			// First set up input stream
			FileInputStream im = new FileInputStream(infile);
			BufferedInputStream bis = new BufferedInputStream(im);
			BufferedOutputStream bos = new BufferedOutputStream(new FileOutputStream(outfile));
		
			// Now compress data
			DeflaterOutputStream compressor = new DeflaterOutputStream(bos);
			
			byte[] buf = new byte[1024];

			int bytesread=0;
			int totalbytesread=0;
			
			while ((bytesread = bis.read(buf))>0){
				totalbytesread=totalbytesread+bytesread;
				compressor.write(buf,0,bytesread);
			}
			compressor.flush();
			compressor.close();

			if(verbose) System.err.println("Read "+totalbytesread+" bytes");
			if(verbose) System.err.println("Output file has length: "+outfile.length()+" bytes");
		}
		catch (Exception e) {
			System.out.println(""+e);
		}
		
	}
}
