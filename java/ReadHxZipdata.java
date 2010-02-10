// Simple java program to read zlib format data as found in 
// Amiramesh files using the HxZip encoding
// (c) Copyright 2010 Gregory Jefferis. All Rights Reserved. 

// Note that this seems to fail silently when decompressing files over AFP
// Really can't think why

// OK more testing reveals silent failures on local disk as well, 
// when decompressing whole buffer at once

// Solution appears to be to read compressed data into memory first


import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.BufferedInputStream;
import java.util.zip.Inflater;

public class ReadHxZipdata {
	
	public static void main(String [ ] args)
	{
		if(args.length!=5) {
			System.err.println("usage: java ReadHxZipdata <infile> <offset> <compressedDatalength> <uncompressedDatalength> <outfile>");
			return;
		}

		String file=args[0];
		String outfile=args[4];
		int offset=0,dataLength=0,compressedDataLength=0;
		try {
			// the String to int conversion happens here
			offset = Integer.parseInt(args[1].trim());
			compressedDataLength = Integer.parseInt(args[2].trim());
			dataLength = Integer.parseInt(args[3].trim());
			System.err.println("dataLength is "+dataLength+" bytes");
	    }
	    catch (NumberFormatException nfe)
	    {
			System.out.println("NumberFormatException: " + nfe.getMessage());
	    }
		
		try{
			
			// First read compressed data
			FileInputStream im = new FileInputStream(file);
			System.err.println("Skipping "+offset+" bytes");
			BufferedInputStream bis = new BufferedInputStream(im);
			bis.skip(offset);
			byte[] compressedData = new byte[compressedDataLength];
			bis.read(compressedData);
			im.close();

			// Now decompress data
			Inflater decompressor = new Inflater();
			byte[] buf = new byte[dataLength];
			decompressor.setInput(compressedData, 0, compressedDataLength);
			int bytesread = decompressor.inflate(buf);
			decompressor.end();
			
			System.err.println("Read "+bytesread+" bytes");
			
			// Now write output
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
