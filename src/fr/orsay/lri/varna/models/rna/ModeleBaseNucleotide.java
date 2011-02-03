/*
 VARNA is a tool for the automated drawing, visualization and annotation of the secondary structure of RNA, designed as a companion software for web servers and databases.
 Copyright (C) 2008  Kevin Darty, Alain Denise and Yann Ponty.
 electronic mail : Yann.Ponty@lri.fr
 paper mail : LRI, bat 490 Université Paris-Sud 91405 Orsay Cedex France

 This file is part of VARNA version 3.1.
 VARNA version 3.1 is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

 VARNA version 3.1 is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along with VARNA version 3.1.
 If not, see http://www.gnu.org/licenses.
 */
package fr.orsay.lri.varna.models.rna;

import java.awt.geom.Point2D;


/**
 * The rna base model with the first character of the nitrogenous base and it
 * display
 * 
 * @author darty
 * 
 */
public class ModeleBaseNucleotide extends ModeleBase {

	/**
	 * 
	 */
	private static final long serialVersionUID = -5493938366569588113L;
	private Character _c;
	private int _index;

	/**
	 * Creates a new rna base with the default display style and a space as
	 * nitrogenous base
	 */
	public ModeleBaseNucleotide(int index) {
		this(' ', index);
	}

	/**
	 * Creates a new rna base with the nitrogenous base
	 * 
	 * @param c
	 *            The first character of the nitrogenous base
	 */
	public ModeleBaseNucleotide(char c, int index) {
		this(c, new ModeleStyleBase(), index);
	}

	/**
	 * Creates a new rna base with the nitrogenous base
	 * 
	 * @param c
	 *            The first character of the nitrogenous base
	 */
	public ModeleBaseNucleotide(char c, int index, int baseNumber) {
		this(c, new ModeleStyleBase(), index);
		_realIndex = baseNumber;
	}
	
	/**
	 * Creates a new rna base with the nitrogenous base and the display style
	 * 
	 * @param c
	 *            The first character of the nitrogenous base
	 * @param msb
	 *            The display style
	 */
	public ModeleBaseNucleotide(char c, ModeleStyleBase msb, int index) {
		this(new Point2D.Double(), new Point2D.Double(), true, c, msb, -1,
				index);
	}

	/**
	 * Creates a new rna base with a display style
	 * 
	 * @param msb
	 *            The display style
	 */
	public ModeleBaseNucleotide(ModeleStyleBase msb, int index, int baseNumber) {
		this(' ', msb, index);
		_realIndex = baseNumber;
	}

	/**
	 * Creates a new rna base with a space as the nitrogenous base and the
	 * display style
	 * 
	 * @param coord
	 * @param index
	 */
	public ModeleBaseNucleotide(Point2D.Double coord, int index) {
		this(new Point2D.Double(coord.getX(), coord.getY()),
				new Point2D.Double(), true, ' ', new ModeleStyleBase(), -1,
				index);
	}

	/**
	 * Creates a new rna base from another one with the same attributes
	 * 
	 * @param mb
	 *            The base to copy
	 */
	public ModeleBaseNucleotide(ModeleBaseNucleotide mb, int index) {
		this(
				new Point2D.Double(mb.getCoords().getX(), mb.getCoords()
						.getY()), new Point2D.Double(mb.getCenter().getX(), mb
						.getCenter().getY()), true, mb.get_c(), mb
						.getStyleBase(), mb.getElementStructure(), index);
	}

	public ModeleBaseNucleotide(Point2D.Double coords, Point2D.Double center,
			boolean colorie, char label, ModeleStyleBase mb, int elementStruct,
			int index) {
		_colorie = colorie;
		_c = label;
		_styleBase = mb;
		_coords = new VARNAPoint(coords);
		_center = new VARNAPoint(center);
		_elementStructure = elementStruct;
		_index = index;
		_realIndex = index + 1;
		_value = 0.0;
	}

	public ModeleStyleBase getStyleBase() {
		if (_colorie)
			return _styleBase;
		return new ModeleStyleBase();
	}

	public Character get_c() {
		return _c;
	}

	public void set_c(Character _c) {
		this._c = _c;
	}

	public Boolean getColorie() {
		return _colorie;
	}

	public void setColorie(Boolean _colorie) {
		this._colorie = _colorie;
	}


	public int getElementStructure() {
		return _elementStructure;
	}

	public String getContent() {
		return "" + _c;
	}

	public int getIndex() {
		return _index;
	}
	
	public String toString()
	{
		return ""+this._realIndex+" ("+_index+") (x,y):"+this._coords +" C:"+_center;
	}

}
