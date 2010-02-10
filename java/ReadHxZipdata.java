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
import java.util.zip.InflaterInputStream;

public class ReadHxZipdata {
	
	public static void main(String [ ] args)
	{
		if(args.length!=4) {
			System.err.println("usage: java ReadHxZipdata <infile> <zlibdataoffset> <uncompressedDatalength> <outfile>");
			return;
		}

		String file=args[0];
		String outfile=args[3];
		int zlibdataoffset=0,dataLength=0,compressedDataLength=0;
		try {
			// the String to int conversion happens here
			zlibdataoffset = Integer.parseInt(args[1].trim());
			dataLength = Integer.parseInt(args[2].trim());
			System.err.println("dataLength is "+dataLength+" bytes");
	    }
	    catch (NumberFormatException nfe)
	    {
			System.out.println("NumberFormatException: " + nfe.getMessage());
	    }
		
		try{
			
			// First set up input stream
			FileInputStream im = new FileInputStream(file);
			System.err.println("Skipping "+zlibdataoffset+" bytes");
			BufferedInputStream bis = new BufferedInputStream(im);
			bis.skip(zlibdataoffset);
		
			// Now decompress data
			InflaterInputStream decompressor = new InflaterInputStream(bis);
			byte[] buf = new byte[dataLength];
			
			// Loop to read uncompressed data until we've got all we wanted
		    int offset = 0, len = dataLength;
			while (len > 0 && decompressor.available()==1) {
				int count = decompressor.read(buf, offset, len);
				if (count <= 0)
					return;
				offset += count;
				len -= count;
			}
			im.close();
						
			int bytesread = dataLength - len;
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
