/*
 VARNA is a tool for the automated drawing, visualization and annotation of the secondary structure of RNA, designed as a companion software for web servers and databases.
 Copyright (C) 2008  Kevin Darty, Alain Denise and Yann Ponty.
 electronic mail : Yann.Ponty@lri.fr
 paper mail : LRI, bat 490 Universit� Paris-Sud 91405 Orsay Cedex France

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

import java.awt.Point;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.StringReader;
import java.net.URL;
import java.net.URLConnection;
import java.util.Stack;
import java.util.Vector;

import org.xml.sax.Attributes;
import org.xml.sax.InputSource;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.DefaultHandler;

public class RNAMLParser extends DefaultHandler {
	public class HelixTemp {
		public int pos5, pos3, length;
		public String name;

		public HelixTemp(int pos5, int pos3, int length, String name) {
			this.pos3 = pos3;
			this.pos5 = pos5;
			this.length = length;
			this.name = name;
		}

		public String toString() {
			return ("[" + name + "," + pos5 + "," + pos3 + "," + length + "]");
		}

	}

	public class BPTemp {
		public int pos5, pos3;
		public String edge5, edge3, orientation;

		public BPTemp(int pos5, int pos3, String edge5, String edge3,
				String orientation) {
			this.pos5 = pos5;
			this.pos3 = pos3;
			this.edge5 = edge5;
			this.edge3 = edge3;
			this.orientation = orientation;
		}

		public ModeleStyleBP createBPStyle(ModeleBase mb5, ModeleBase mb3) {
			ModeleStyleBP.Edge e5, e3;
			@SuppressWarnings("unused")
			boolean isCanonical = false;
			if (edge5.equals("W")) {
				e5 = ModeleStyleBP.Edge.WATSON_CRICK;
			} else if (edge5.equals("H")) {
				e5 = ModeleStyleBP.Edge.HOOGSTEEN;
			} else if (edge5.equals("S")) {
				e5 = ModeleStyleBP.Edge.SUGAR;
			} else {
				e5 = ModeleStyleBP.Edge.WATSON_CRICK;
			}

			if (edge3.equals("W")) {
				e3 = ModeleStyleBP.Edge.WATSON_CRICK;
			} else if (edge3.equals("H")) {
				e3 = ModeleStyleBP.Edge.HOOGSTEEN;
			} else if (edge3.equals("S")) {
				e3 = ModeleStyleBP.Edge.SUGAR;
			} else {
				e3 = ModeleStyleBP.Edge.WATSON_CRICK;
			}

			if ((edge5.equals("+") && edge3.equals("+"))
					|| (edge5.equals("-") && edge3.equals("-"))) {
				e3 = ModeleStyleBP.Edge.WATSON_CRICK;
				e5 = ModeleStyleBP.Edge.WATSON_CRICK;
			}

			ModeleStyleBP.Stericity ster;

			if (orientation.equals("c")) {
				ster = ModeleStyleBP.Stericity.CIS;
			} else if (orientation.equals("t")) {
				ster = ModeleStyleBP.Stericity.TRANS;
			} else {
				ster = ModeleStyleBP.Stericity.CIS;
			}

			return (new ModeleStyleBP(mb5, mb3, e5, e3, ster));
		}

		public String toString() {
			return ("[" + pos5 + "," + pos3 + "," + edge5 + "," + edge3 + ","
					+ orientation + "]");
		}
	}

	private String _sequence = "";
	private Vector<Integer> _sequenceIDs = new Vector<Integer>();
	private Vector<BPTemp> _structure = new Vector<BPTemp>();
	private Vector<HelixTemp> _helices = new Vector<HelixTemp>();

	@SuppressWarnings("unused")
	private boolean _inSequenceIDs, _inLength, _inSequence, _inHelix,
			_inStrAnnotation, _inBP, _inBP5, _inBP3, _inEdge5, _inEdge3,
			_inPosition, _inBondOrientation;
	private StringBuffer _buffer;
	private String _firstModel = "";
	private String _currentModel = "";
	private int _id5, _id3, _length;
	String _edge5, _edge3, _orientation, _helixID;

	private Vector<BPTemp> _structurePlanar = new Vector<BPTemp>();
	private Vector<BPTemp> _structureAux = new Vector<BPTemp>();

	public RNAMLParser() {
		super();
		_inSequenceIDs = false;
		_inSequence = false;
		_inStrAnnotation = false;
		_inBP = false;
		_inBP5 = false;
		_inBP3 = false;
		_inPosition = false;
		_inEdge5 = false;
		_inEdge3 = false;
		_inBondOrientation = false;
		_inHelix = false;
	}

	public InputSource createSourceFromURL(String path)
	{
		URL url = null;
		try {
			url = new URL(path);
			URLConnection connexion = url.openConnection();
			connexion.setUseCaches(false);
			InputStream r = connexion.getInputStream();
			InputStreamReader inr = new InputStreamReader(r);
			return new InputSource(inr);
		}
		catch(Exception e)
		{
			e.printStackTrace();
		}
		return new InputSource(new StringReader(""));
	}
	
	public InputSource resolveEntity(String publicId, String systemId) {
		System.out.println("[crade]");
		if (systemId.endsWith("rnaml.dtd"))
			return createSourceFromURL("http://varna.lri.fr/bin/rnaml.dtd");
		else
			return new InputSource(new StringReader(""));
	}

	public void startElement(String uri, String localName, String qName,
			Attributes attributes) throws SAXException {
		if (qName.equals("numbering-table")) {
			_inSequenceIDs = true;
			_buffer = new StringBuffer();
		} else if (qName.equals("helix")) {
			_inHelix = true;
			_buffer = new StringBuffer();
			_helixID = attributes.getValue("id");
		} else if (qName.equals("seq-data")) {
			_inSequence = true;
			_buffer = new StringBuffer();
		} else if (qName.equals("length")) {
			_inLength = true;
			_buffer = new StringBuffer();
		} else if (qName.equals("str-annotation")) {
			_inStrAnnotation = true;
		} else if (qName.equals("base-pair")) {
			_inBP = true;
		} else if (qName.equals("base-id-5p")) {
			if (_inBP || _inHelix) {
				_inBP5 = true;
			}
		} else if (qName.equals("base-id-3p")) {
			if (_inBP || _inHelix) {
				_inBP3 = true;
			}
		} else if (qName.equals("edge-5p")) {
			_inEdge5 = true;
			_buffer = new StringBuffer();
		} else if (qName.equals("edge-3p")) {
			_inEdge3 = true;
			_buffer = new StringBuffer();
		} else if (qName.equals("position")) {
			_inPosition = true;
			_buffer = new StringBuffer();
		} else if (qName.equals("bond-orientation")) {
			_inBondOrientation = true;
			_buffer = new StringBuffer();
		} else if (qName.equals("molecule")) {
			String id = (attributes.getValue("id"));
			if (_firstModel.length()==0) {
				_firstModel = id;
			}
			_currentModel = id;
		} else {
			// We don't care too much about the rest ...
		}
	}

	public void endElement(String uri, String localName, String qName)
			throws SAXException {
		if (qName.equals("numbering-table")) {
			_inSequenceIDs = false;
			String content = _buffer.toString();
			content = content.trim();
			String[] tokens = content.split("\\s+");
			Vector<Integer> results = new Vector<Integer>();
			for (int i = 0; i < tokens.length; i++) {
				try {
					results.add(new Integer(Integer.parseInt(tokens[i])));
				} catch (NumberFormatException e) {
					e.printStackTrace();
				}
			}
			if (_currentModel == _firstModel) {
				_sequenceIDs = results;
			}
			_buffer = null;
		} else if (qName.equals("seq-data")) {
			_inSequence = false;
			String content = _buffer.toString();
			content = content.trim();
			String[] tokens = content.split("\\s+");
			StringBuffer results = new StringBuffer();
			for (int i = 0; i < tokens.length; i++) {
				results.append(tokens[i]);
			}

			if (_currentModel == _firstModel) {
				_sequence = results.toString();
			}
			_buffer = null;
		} else if (qName.equals("bond-orientation")) {
			_inBondOrientation = false;
			String content = _buffer.toString();
			content = content.trim();
			_orientation = content;
			_buffer = null;
		} else if (qName.equals("str-annotation")) {
			_inStrAnnotation = false;
		} else if (qName.equals("base-pair")) {
			_inBP = false;
			if (_currentModel == _firstModel) {
				BPTemp bp = new BPTemp(_id5, _id3, _edge5, _edge3,_orientation);
				this._structure.add(bp);
			}
		} else if (qName.equals("helix")) {
			_inHelix = false;
			HelixTemp h = new HelixTemp(_id5, _id3, _length, _helixID);
			_helices.add(h);
		} else if (qName.equals("base-id-5p")) {
			_inBP5 = false;
		} else if (qName.equals("base-id-3p")) {
			_inBP3 = false;
		} else if (qName.equals("length")) {
			_inLength = false;
			String content = _buffer.toString();
			content = content.trim();
			_length = Integer.parseInt(content);
			_buffer = null;
		} else if (qName.equals("position")) {
			String content = _buffer.toString();
			content = content.trim();
			int pos = Integer.parseInt(content);
			if (_inBP5) {
				_id5 = pos;
			}
			if (_inBP3) {
				_id3 = pos;
			}
			_buffer = null;
		} else if (qName.equals("edge-5p")) {
			_inEdge5 = false;
			String content = _buffer.toString();
			content = content.trim();
			_edge5 = content;
			_buffer = null;
		} else if (qName.equals("edge-3p")) {
			_inEdge3 = false;
			String content = _buffer.toString();
			content = content.trim();
			_edge3 = content;
			_buffer = null;
		} else {
			// We don't care too much about the rest ...
		}
	}

	public void characters(char[] ch, int start, int length)
			throws SAXException {
		String lecture = new String(ch, start, length);
		if (_buffer != null)
			_buffer.append(lecture);
	}

	public void startDocument() throws SAXException {
	}

	public void endDocument() throws SAXException {
		postProcess();
	}

	private void filterBasePairs() {
		Vector<BPTemp> result = new Vector<BPTemp>();
		for (int i = 0; i < _structure.size(); i++) {
			BPTemp bp = _structure.get(i);
			if (bp.orientation.equals("c") || bp.orientation.equals("t")) {
				result.add(bp);
			}
		}
		_structure = result;
	}

	public String getSequence() {
		return _sequence;
	}

	public static boolean isSelfCrossing(int[] str) {
		Stack<Point> intervals = new Stack<Point>();
		intervals.add(new Point(0, str.length - 1));
		while (!intervals.empty()) {
			Point p = intervals.pop();
			if (p.x <= p.y) {
				if (str[p.x] == -1) {
					intervals.push(new Point(p.x + 1, p.y));
				} else {
					int i = p.x;
					int j = p.y;
					int k = str[i];
					if ((k <= i) || (k > j)) {
						return true;
					} else {
						intervals.push(new Point(i + 1, k - 1));
						intervals.push(new Point(k + 1, j));
					}
				}
			}
		}
		return false;
	}

	@SuppressWarnings("unused")
	private void debugPrintArray(int[] str) {
		StringBuffer s = new StringBuffer("[");
		for (int i = 0; i < str.length; i++) {
			if (i != 0) {
				s.append(",");
			}
			s.append(str[i]);

		}
		s.append("]");
		System.out.println(s.toString());
	}

	/**
	 * Computes and returns a maximal planar subset of the current structure. 
	 * @param str A sequence of base-pairing positions
	 * @return A sequence of non-crossing base-pairing positions
	 */
	
	public static int[] planarize(int[] str) {
		if (!isSelfCrossing(str)) {
			return str;
		}
		
		int length = str.length;

		int[] result = new int[length];
		for (int i = 0; i < result.length; i++) {
			result[i] = -1;
		}

		short[][] tab = new short[length][length];
		short[][] backtrack = new short[length][length];
		int theta = 3;

		for (int i = 0; i < result.length; i++) {
			for (int j = i; j < Math.min(i + theta, result.length); j++) {
				tab[i][j] = 0;
				backtrack[i][j] = -1;
			}
		}
		for (int n = theta; n < length; n++) {
			for (int i = 0; i < length - n; i++) {
				int j = i + n;
				tab[i][j] = tab[i + 1][j];
				backtrack[i][j] = -1;
				int k = str[i];
				if ((k != -1) && (k <= j) && (i < k)) {
					int tmp = 1;
					if (i + 1 <= k - 1) {
						tmp += tab[i + 1][k - 1];
					}
					if (k + 1 <= j) {
						tmp += tab[k + 1][j];
					}
					if (tmp > tab[i][j]) {
						tab[i][j] = (short) tmp;
						backtrack[i][j] = (short) k;
					}
				}
			}
		}
		Stack<Point> intervals = new Stack<Point>();
		intervals.add(new Point(0, length - 1));
		while (!intervals.empty()) {
			Point p = intervals.pop();
			if (p.x <= p.y) {
				if (backtrack[p.x][p.y] == -1) {
					result[p.x] = -1;
					intervals.push(new Point(p.x + 1, p.y));
				} else {
					int i = p.x;
					int j = p.y;
					int k = backtrack[p.x][p.y];
					result[i] = k;
					result[k] = i;
					intervals.push(new Point(i + 1, k - 1));
					intervals.push(new Point(k + 1, j));
				}
			}
		}
		return result;
	}

	public int[] getBasicPlanarStructure() {
		int[] str = new int[_sequence.length()];
		for (int i = 0; i < str.length; i++) {
			str[i] = -1;
		}
		// System.out.print("Seq:"+_sequence+"\nHel:"+_helices);

        // First for the canonical bps
		for (int i = 0; i < _structure.size(); i++) {
			BPTemp bp = _structure.get(i);
			int a = bp.pos5 - 1;
			int b = bp.pos3 - 1;
			if ((str[a]==-1)&&(str[b]==-1)
					&& (bp.orientation.equals("c"))
					&&(    (bp.edge5.equals("+") && bp.edge3.equals("+"))
					    || (bp.edge5.equals("-") && bp.edge3.equals("-")))) {
				str[a] = b;
				str[b] = a;
			}
		}
		
		// ... then helices...
		 for (int i=0;i<_helices.size();i++) { 
			 HelixTemp h = _helices.get(i);
			 for (int j=0;j<h.length;j++) { 
				 int a = h.pos5-1+j; 
				 int b = h.pos3-1-j; 
				 if ((str[a]==-1)&&(str[b]==-1))
				 {
				 str[a] = b; 
				 str[b] = a; 
				}
			 } 
		}
		
		// ... and finally non-canonical bps
		for (int i = 0; i < _structure.size(); i++) {
			Vector<Integer> basenumbers=getBaseNumbers();
			BPTemp bp = _structure.get(i);
			int a = bp.pos5 - 1;
			int b = bp.pos3 - 1;
			
			int realPosA= basenumbers.get(a);
			int realPosB= basenumbers.get(b);
			//System.out.print("a:"+a+" real:"+realPosA+"\n");
			if ((str[a]==-1)&&(str[b]==-1)&&(realPosB-realPosA>4)) {
				str[a] = b;
				str[b] = a;
			}
		}
		/*for(int j=0;j<str.length;j++){
			System.out.print(str[j]+"  ");
		}*/
		//System.out.print("\n");
		return planarize(str);
	}

	@SuppressWarnings("unused")
	private BPTemp[] getBPArray() {
		BPTemp[] result = new BPTemp[_sequence.length()];
		for (int i = 0; i < result.length; i++) {
			result[i] = null;
		}

		Vector<BPTemp> res = _structure;
		for (int i = 0; i < res.size(); i++) {
			BPTemp bp = res.get(i);
			int bpfrom = bp.pos5 - 1;
			int bpto = bp.pos3 - 1;
			if ((result[bpfrom] == null) && (result[bpto] == null)) {
				result[bpfrom] = bp;
				result[bpto] = bp;
			}
		}
		return result;
	}

	private void postProcess() {
		// First, check if base numbers were specified
		if (_sequenceIDs.size()==0)
		{
			Vector<Integer> results = new Vector<Integer>();
			for (int i = 0; i < _sequence.length(); i++) {
				results.add(new Integer(i+1));
			}
			_sequenceIDs = results;
		}
		filterBasePairs();
		int[] planar = getBasicPlanarStructure();
		for (int i = 0; i < _structure.size(); i++) {
			BPTemp bp = _structure.get(i);
			int k = bp.pos5 - 1;
			int l = bp.pos3 - 1;
			if ((k == i) && (l == planar[i])) {
				_structurePlanar.add(bp);
			} else {
				_structureAux.add(bp);
			}
		}
	}

	public Vector<BPTemp> getPlanarBPs() {
		return _structurePlanar;
	}

	public Vector<BPTemp> getAuxBPs() {
		return _structureAux;
	}

	public Vector<Integer> getBaseNumbers() {
		return _sequenceIDs;
	}

}