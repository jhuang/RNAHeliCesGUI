package fr.orsay.lri.varna.models;

import java.awt.Color;
import java.util.ArrayList;
import java.util.Collection;

import fr.orsay.lri.varna.models.rna.ModeleBase;

public class BaseList {
	private ArrayList<ModeleBase> _bases = new ArrayList<ModeleBase>(); 
	private String _caption;

	public BaseList( String caption)
	{
		_caption = caption;
	}
	
	
	public BaseList( String caption, ModeleBase mb)
	{
		_caption = caption;
		addBase(mb);
	}

	public boolean contains(ModeleBase mb)
	{
		return _bases.contains(mb);
	}

	
	public String getCaption()
	{
		return _caption;
	}
	
	public void addBase(ModeleBase b)
	{
		_bases.add(b);
	}

	public void removeBase(ModeleBase b)
	{
		_bases.remove(b);
	}

	
	public void addBases(Collection<? extends ModeleBase> mbs)
	{
		_bases.addAll(mbs);
	}

	public ArrayList<ModeleBase> getBases()
	{
		return _bases;
	}

	public void clear()
	{
		_bases.clear();
	}

	public static Color getAverageColor(ArrayList<Color> cols)
	{
		int r=0,g=0,b=0;
		for (Color c : cols)
		{
			r += c.getRed();
			g += c.getGreen();
			b += c.getBlue();
		}
		if (cols.size()>0)
		{ 
			r /= cols.size();
			g /= cols.size();
			b /= cols.size();
		}
		return new Color(r,g,b);
	}
	
	public Color getAverageOutlineColor()
	{
		ArrayList<Color> cols = new ArrayList<Color>(); 
		for (ModeleBase mb : _bases)
		{  cols.add(mb.getStyleBase().get_base_outline_color()); }
		return getAverageColor(cols);
	}

	public Color getAverageNameColor()
	{
		ArrayList<Color> cols = new ArrayList<Color>(); 
		for (ModeleBase mb : _bases)
		{  cols.add(mb.getStyleBase().get_base_name_color()); }
		return getAverageColor(cols);
	}

	public Color getAverageNumberColor()
	{
		ArrayList<Color> cols = new ArrayList<Color>(); 
		for (ModeleBase mb : _bases)
		{  cols.add(mb.getStyleBase().get_base_number_color()); }
		return getAverageColor(cols);
	}

	public Color getAverageInnerColor()
	{
		ArrayList<Color> cols = new ArrayList<Color>(); 
		for (ModeleBase mb : _bases)
		{  cols.add(mb.getStyleBase().get_base_inner_color()); }
		return getAverageColor(cols);
	}

	public String getNumbers()
	{
		String result = ""; 
		for (int i=0; i<_bases.size();i++)
		{  
			if (i>0)
			{result += ",";}
			ModeleBase mb  = _bases.get(i);
			result += "" + mb.getBaseNumber(); 
		}
		result += "";
		return result;
	}

	public String getContents()
	{
		String result = ""; 
		for (int i=0; i<_bases.size();i++)
		{  
			if (i>0)
			{result += ",";}
			ModeleBase mb  = _bases.get(i);
			result += mb.getContent(); 
		}
		result += "";
		return result;
	}

	public ArrayList<Integer> getIndices()
	{
		ArrayList<Integer> indices = new ArrayList<Integer>();
		for (ModeleBase mb : _bases)
		{
			indices.add(mb.getIndex());
		}
		return indices;
	}
}
