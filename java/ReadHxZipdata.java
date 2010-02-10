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
			System.err.println("dataLength is "+dataLength+" bytes");
	    }
	    catch (NumberFormatException nfe)
	    {
			System.out.println("NumberFormatException: " + nfe.getMessage());
	    }
		
		try{
			
			FileInputStream im = new FileInputStream(file);
			im.skip(offset);
			System.err.println("Skipping "+offset+" bytes");
			InflaterInputStream decompressor = new InflaterInputStream(im);
			byte[] buf = new byte[dataLength];
			int bytesread = decompressor.read(buf,0,dataLength);
			im.close();
			decompressor.close();
			
			System.err.println("Read "+bytesread+" bytes");
			
			FileOutputStream fos = new FileOutputStream(outfile);
			fos.write(buf);
			fos.close();
			System.err.println("Wrote "+bytesread+" bytes");
		}
		catch (Exception e) {
			System.out.println(""+e);
		}
		
	}
}
