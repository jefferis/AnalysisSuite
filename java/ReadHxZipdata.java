// Simple java program to read zlib format data as found in 
// Amiramesh files using the HxZip encoding
// (c) Copyright 2010 Gregory Jefferis. All Rights Reserved. 

import java.io.UnsupportedEncodingException;
import java.io.PrintWriter;
import java.io.OutputStreamWriter;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.util.zip.InflaterInputStream;

public class ReadHxZipdata {
	private static PrintWriter out;
	private static String charsetName = "UTF-8";
	
	static {
        try {
            out = new PrintWriter(new OutputStreamWriter(System.out, charsetName), true);
        }
        catch (UnsupportedEncodingException e) { System.out.println(e); }
    }

	public static void main(String [ ] args)
	{
		if(args.length!=4) {
			System.err.println("usage: java ReadHxZipdata <infile> <offset> <uncompressedDatalength> <outfile>");
			return;
		}

		String file=args[0];
		String outfile=args[3];
		int offset=0,dataLength=0;
		try {
			// the String to int conversion happens here
			offset = Integer.parseInt(args[1].trim());
			dataLength = Integer.parseInt(args[2].trim());
	    }
	    catch (NumberFormatException nfe)
	    {
			System.out.println("NumberFormatException: " + nfe.getMessage());
	    }
		
		try{
			FileOutputStream fos = new FileOutputStream(outfile);
			
			FileInputStream im = new FileInputStream(file);
			im.skip(offset);
			InflaterInputStream decompressor = new InflaterInputStream(im);
			byte[] buf = new byte[dataLength];
			decompressor.read(buf,0,dataLength);
			
			fos.write(buf);
		}
		catch (Exception e) {
			System.out.println(""+e);
		}
		
	}
}
