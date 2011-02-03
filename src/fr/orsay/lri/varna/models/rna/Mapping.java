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

import java.io.Serializable;
import java.util.Enumeration;
import java.util.Hashtable;

import fr.orsay.lri.varna.exceptions.MappingException;

public class Mapping implements Serializable {

	/**
	 * 
	 */
	private static final long serialVersionUID = -3031358968555310380L;

	public static final int UNKNOWN = -1;

	Hashtable<Integer, Integer> _mapping = new Hashtable<Integer, Integer>();
	Hashtable<Integer, Integer> _invMapping = new Hashtable<Integer, Integer>();

	public Mapping() {

	}

	public void addCouple(int i, int j) throws MappingException {
		if (_mapping.containsKey(i) || _invMapping.containsKey(j)) {
			throw new MappingException(
					MappingException.MULTIPLE_PARTNERS_DEFINITION_ATTEMPT);
		}

		_mapping.put(new Integer(i), new Integer(j));
		_invMapping.put(new Integer(j), new Integer(i));
	}

	public int getPartner(int i) {
		if (!_mapping.containsKey(i))
			return UNKNOWN;
		else
			return _mapping.get(i);
	}

	public int getAncestor(int j) {
		if (!_invMapping.containsKey(j))
			return UNKNOWN;
		else
			return _invMapping.get(j);
	}

	public int[] getSourceElems() {
		int[] elems = new int[_mapping.size()];
		Enumeration<Integer> en = _mapping.keys();
		int i = 0;
		while (en.hasMoreElements()) {
			int a = en.nextElement();
			elems[i] = a;
			i++;
		}
		return elems;
	}

	public int[] getTargetElems() {
		int[] elems = new int[_invMapping.size()];
		Enumeration<Integer> en = _invMapping.keys();
		int i = 0;
		while (en.hasMoreElements()) {
			int a = en.nextElement();
			elems[i] = a;
			i++;
		}
		return elems;
	}

	public static Mapping readMappingFromAlignment(String m, String n)
			throws MappingException {
		Mapping map = new Mapping();
		if (m.length() != n.length()) {
			throw new MappingException(MappingException.BAD_ALIGNMENT_INPUT);
		}
		int i = 0;
		int j = 0;
		for (int k = 0; k < m.length(); k++) {
			char a = m.charAt(k);
			char b = n.charAt(k);
			if ((a != '-') && (a != ':') && (b != '-') && (b != ':')) {
				map.addCouple(i, j);
			}

			if ((a != '-') && (a != ':')) {
				j++;
			}
			if ((b != '-') && (b != ':')) {
				i++;
			}
		}
		return map;
	}

	public static Mapping DefaultMapping(int n, int m) {
		Mapping map = new Mapping();
		try {
			for (int i = 0; i < Math.min(n, m); i++) {
				map.addCouple(i, i);
			}
		} catch (MappingException e) {
			e.printStackTrace();
		}
		return map;
	}

	public static Mapping DefaultOutermostMapping(int n, int m) {
		Mapping map = new Mapping();
		try {
			int k = Math.min(n, m);
			int i = 0;
			int j = 0;
			boolean pile = true;
			while (i <= (k - 1) - j) {
				if (pile) {
					map.addCouple(i, i);
					i++;
				} else {
					map.addCouple(n - 1 - j, m - 1 - j);
					j++;
				}
				pile = !pile;
			}
		} catch (MappingException e) {
			e.printStackTrace();
		}
		// System.println(map);
		return map;
	}

	public String toString() {
		String tmp = "";
		Enumeration<Integer> en = _mapping.keys();
		while (en.hasMoreElements()) {
			Integer i = en.nextElement();
			Integer j = _mapping.get(i);
			tmp += "(" + i + "," + j + ") ";
		}
		return tmp;
	}
}
