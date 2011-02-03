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

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Point;
import java.awt.Shape;
import java.awt.Stroke;
import java.awt.geom.GeneralPath;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Reader;
import java.io.Serializable;
import java.io.StreamTokenizer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.Stack;
import java.util.Vector;

import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;

import org.xml.sax.InputSource;

import fr.orsay.lri.varna.exceptions.ExceptionExportFailed;
import fr.orsay.lri.varna.exceptions.ExceptionFileFormatOrSyntax;
import fr.orsay.lri.varna.exceptions.ExceptionInvalidRNATemplate;
import fr.orsay.lri.varna.exceptions.ExceptionLoadingFailed;
import fr.orsay.lri.varna.exceptions.ExceptionNAViewAlgorithm;
import fr.orsay.lri.varna.exceptions.ExceptionPermissionDenied;
import fr.orsay.lri.varna.exceptions.ExceptionUnmatchedClosingParentheses;
import fr.orsay.lri.varna.exceptions.ExceptionWritingForbidden;
import fr.orsay.lri.varna.interfaces.InterfaceVARNAListener;
import fr.orsay.lri.varna.interfaces.InterfaceVARNAObservable;
import fr.orsay.lri.varna.models.CubicBezierCurve;
import fr.orsay.lri.varna.models.VARNAConfig;
import fr.orsay.lri.varna.models.annotations.ChemProbAnnotation;
import fr.orsay.lri.varna.models.annotations.HighlightRegionAnnotation;
import fr.orsay.lri.varna.models.annotations.TextAnnotation;
import fr.orsay.lri.varna.models.export.PSExport;
import fr.orsay.lri.varna.models.export.SVGExport;
import fr.orsay.lri.varna.models.export.SecStrDrawingProducer;
import fr.orsay.lri.varna.models.export.XFIGExport;
import fr.orsay.lri.varna.models.naView.NAView;
import fr.orsay.lri.varna.models.templates.RNATemplate;
import fr.orsay.lri.varna.models.templates.RNATemplateAlign;
import fr.orsay.lri.varna.models.templates.RNATemplateDrawingAlgorithmException;
import fr.orsay.lri.varna.models.templates.RNATemplateMapping;
import fr.orsay.lri.varna.models.templates.RNATemplate.EdgeEndPointPosition;
import fr.orsay.lri.varna.models.templates.RNATemplate.In1Is;
import fr.orsay.lri.varna.models.templates.RNATemplate.RNATemplateElement;
import fr.orsay.lri.varna.models.templates.RNATemplate.RNATemplateHelix;
import fr.orsay.lri.varna.models.templates.RNATemplate.RNATemplateUnpairedSequence;
import fr.orsay.lri.varna.models.templates.RNATemplate.RNATemplateElement.EdgeEndPoint;

import java.lang.Math;


/**
 * The RNA model which contain the base list and the draw algorithm mode
 * 
 * @author darty
 * 
 */
public class RNA extends InterfaceVARNAObservable implements Serializable{
	/**
	 * 
	 */
	private static final long serialVersionUID = 7541274455751497303L;
	/**
	 * Selects the "Feynman diagram" drawing algorithm that places the bases on
	 * a circle and draws the base-pairings as chords of the circle graph.
	 */

	public static final int FILE_TYPE_BPSEQ = 1;
	public static final int FILE_TYPE_CT = 2;
	public static final int FILE_TYPE_DBN = 3;
	public static final int FILE_TYPE_RNAML = 4;
	public static final int FILE_TYPE_UNKNOWN = 5;
	
	
	
	
	public static final int DRAW_MODE_CIRCULAR = 1;
	/**
	 * Selects the "tree drawing" algorithm. Draws each loop on a circle whose
	 * radius depends on the number of bases involved in the loop. As some
	 * helices can be overlapping in the result, basic interaction is provided
	 * so that the user can "disentangle" the drawing by spinning the helices
	 * around the axis defined by their multiloop (bulge or internal loop)
	 * origin. This is roughly the initial placement strategy of RNAViz.
	 * 
	 * @see <a href="http://rnaviz.sourceforge.net/">RNAViz</a>
	 */
	public static final int DRAW_MODE_RADIATE = 2;
	
	/**
	 * Selects the NAView algorithm.
	 */
	public static final int DRAW_MODE_NAVIEW = 3;
	/**
	 * Selects the linear algorithm.
	 */
	public static final int DRAW_MODE_LINEAR = 4;

	public static final int DRAW_MODE_VARNA_VIEW = 5;
	
	/**
	 * Selects the RNAView algorithm.
	 */
	public static final int DRAW_MODE_MOTIFVIEW = 6;
	
	public static final int DRAW_MODE_TEMPLATE = 7;

	public static final int DEFAULT_DRAW_MODE = DRAW_MODE_RADIATE;
	
	private Double _spaceBetweenBases = 1.0;
	/**
	 * The draw algorithm mode
	 */
	private int _drawMode = DRAW_MODE_RADIATE;

	public int BASE_RADIUS = 10;
	public double LOOP_DISTANCE = 40.0;
	public double BASE_PAIR_DISTANCE = 65.0;
	public double MULTILOOP_DISTANCE = 35.0;

	public double CHEM_PROB_DIST = 14;
	public double CHEM_PROB_BASE_LENGTH = 30;
	public double CHEM_PROB_ARROW_HEIGHT = 10;
	public double CHEM_PROB_ARROW_WIDTH = 5;
	public double CHEM_PROB_TRIANGLE_WIDTH = 2.5;
	public double CHEM_PROB_PIN_SEMIDIAG = 6;
	public double CHEM_PROB_DOT_RADIUS = 6.;
	
	public GeneralPath _debugShape = null;
	
	
	private boolean _comparisonMode;
	private boolean _drawn = false;
	public double _bpHeightIncrement = VARNAConfig.DEFAULT_BP_INCREMENT;

	/**
	 * the base list
	 */
	protected ArrayList<ModeleBase> _listeBases;
	
	/**
	 * the strand list
	 */
	StructureTemp _listStrands = new StructureTemp();
	
	/**
	 * Additional bonds and info can be specified here.
	 */
	private ArrayList<ModeleStyleBP> _structureAux = new ArrayList<ModeleStyleBP>();

	transient private ArrayList<InterfaceVARNAListener> _listeVARNAListener = new ArrayList<InterfaceVARNAListener>();

	ArrayList<Character> _normalBases = new ArrayList<Character>();
	
	private ArrayList<TextAnnotation> _listeAnnotations = new ArrayList<TextAnnotation>();
	private ArrayList<HighlightRegionAnnotation> _listeRegionHighlights = new ArrayList<HighlightRegionAnnotation>();

	boolean _flatExteriorLoop = false;
	
	private String _name = "";

	public RNA() {
		this("");
	}
	public RNA(String n) {
		_name = n;
		_comparisonMode = false;
		_listeBases = new ArrayList<ModeleBase>();
		_drawn = false;
		init();
	}

	public String toString()
	{
		if (_name.equals(""))
		{
			return getStructDBN();
		}
		else
		{
			return _name; 
		}
	}
	
	public RNA(RNA r) {
		_spaceBetweenBases = r._spaceBetweenBases;
		_drawMode = r._drawMode;
		_comparisonMode = r._comparisonMode;
		_listeBases.addAll(r._listeBases);
		_listeVARNAListener = (ArrayList<InterfaceVARNAListener>) r._listeVARNAListener;
		_drawn = r._drawn;
		init();
	}

	public void init() {
		_normalBases.add('a');
		_normalBases.add('c');
		_normalBases.add('g');
		_normalBases.add('u');
		_normalBases.add('-');
	}

	public void saveRNADBN(String path, String title)
			throws ExceptionWritingForbidden {
		try {
			FileWriter out = new FileWriter(path);
			if (!title.equals("")) {
				out.write("> " + title + "\n");
			}
			out.write(getListeBasesToString());
			out.write('\n');
			String str = "";
			for (int i = 0; i < _listeBases.size(); i++) {
				if (_listeBases.get(i).getElementStructure() == -1) {
					str += '.';
				} else {
					if (_listeBases.get(i).getElementStructure() > i) {
						str += '(';
					} else {
						str += ')';
					}
				}
			}
			out.write(str);
			out.write('\n');
			out.close();
		} catch (IOException e) {
			throw new ExceptionWritingForbidden(e.getMessage());
		}
	}

	public Color getBaseInnerColor(int i, VARNAConfig conf) {
		Color result = _listeBases.get(i).getStyleBase()
				.get_base_inner_color();
		char res = _listeBases.get(i).getContent().charAt(0);
		if (conf._drawColorMap)
		{
			result = conf._cm.getColorForValue(_listeBases.get(i).getValue());			
		}
		else if ((conf._colorSpecialBases && !_normalBases.contains(Character.toLowerCase(res)))) {
			result = conf._specialBasesColor;
		} else if ((conf._colorDashBases && (Character.toLowerCase(res) == '-'))) {
			result = conf._dashBasesColor;
		}
		return result;
	}

	public Color getBaseOuterColor(int i, VARNAConfig conf) {
		Color result = _listeBases.get(i).getStyleBase()
				.get_base_outline_color();
		return result;
	}

	public Color getBaseNameColor(int i, VARNAConfig conf) {
		Color result = _listeBases.get(i).getStyleBase().get_base_name_color();
		return result;
	}

	public Color getBasePairColor(ModeleStyleBP bp, VARNAConfig conf) {
		Color bondColor = conf._bondColor;
		if (conf._useBaseColorsForBPs) {
			bondColor = _listeBases.get(bp.getPartner5().getIndex())
					.getStyleBase().get_base_inner_color();
		}
		if (bp != null) {
			bondColor = bp.getColor(bondColor);
		}
		return bondColor;
	}

	public double getBasePairThickness(ModeleStyleBP bp, VARNAConfig conf) {
		double thickness = bp.getThickness(conf._bpThickness);
		return thickness;
	}

	private void drawSymbol(SecStrDrawingProducer out, double posx,
			double posy, double normx, double normy, double radius,
			boolean isCIS, ModeleStyleBP.Edge e, double thickness) {
		Color bck = out.getCurrentColor();
		switch (e) {
		case WATSON_CRICK:
			if (isCIS) {
				out.fillCircle(posx, posy, (radius / 2.0), thickness, bck);
			} else {
				out.drawCircle(posx, posy, (radius / 2.0), thickness);
				out.fillCircle(posx, posy, (radius / 2.0), thickness,
						Color.white);
			}
			break;
		case HOOGSTEEN: {
			double xtab[] = new double[4];
			double ytab[] = new double[4];
			xtab[0] = posx - radius * normx / 2.0 - radius * normy / 2.0;
			ytab[0] = posy - radius * normy / 2.0 + radius * normx / 2.0;
			xtab[1] = posx + radius * normx / 2.0 - radius * normy / 2.0;
			ytab[1] = posy + radius * normy / 2.0 + radius * normx / 2.0;
			xtab[2] = posx + radius * normx / 2.0 + radius * normy / 2.0;
			ytab[2] = posy + radius * normy / 2.0 - radius * normx / 2.0;
			xtab[3] = posx - radius * normx / 2.0 + radius * normy / 2.0;
			ytab[3] = posy - radius * normy / 2.0 - radius * normx / 2.0;
			if (isCIS) {
				out.fillPolygon(xtab, ytab, bck);
			} else {
				out.drawPolygon(xtab, ytab, thickness);
				out.fillPolygon(xtab, ytab,  Color.white);
			}
		}
			break;
		case SUGAR: {
			double ix = radius * normx / 2.0;
			double iy = radius * normy / 2.0;
			double jx = radius * normy / 2.0;
			double jy = -radius * normx / 2.0;
			double xtab[] = new double[3];
			double ytab[] = new double[3];
			xtab[0] = posx - ix + jx;
			ytab[0] = posy - iy + jy;
			xtab[1] = posx + ix + jx;
			ytab[1] = posy + iy + jy;
			xtab[2] = posx - jx;
			ytab[2] = posy - jy;

			if (isCIS) {
				out.fillPolygon(xtab, ytab, bck);
			} else {
				out.drawPolygon(xtab, ytab, thickness);
				out.fillPolygon(xtab, ytab,  Color.white);
			}
		}
			break;
		}
	}

	private void drawBasePair(SecStrDrawingProducer out, Point2D.Double orig,
			Point2D.Double dest, ModeleStyleBP style, VARNAConfig conf) {
		double dx = dest.x - orig.x;
		double dy = dest.y - orig.y;
		double dist = Math.sqrt((dest.x - orig.x) * (dest.x - orig.x)
				+ (dest.y - orig.y) * (dest.y - orig.y));
		dx /= dist;
		dy /= dist;
		double nx = -dy;
		double ny = dx;
		orig = new Point2D.Double(orig.x+BASE_RADIUS*dx,orig.y+BASE_RADIUS*dy);
		dest = new Point2D.Double(dest.x-BASE_RADIUS*dx,dest.y-BASE_RADIUS*dy);
		if (conf._mainBPStyle == VARNAConfig.BP_STYLE.BP_STYLE_LW) {
			double thickness = getBasePairThickness(style, conf);
			double radiusCircle = ((BASE_PAIR_DISTANCE - BASE_RADIUS) / 5.0);

			if (style.isCanonical()) {
				if (style.isCanonicalGC()) {
					if ((orig.x != dest.x) || (orig.y != dest.y)) {
						nx *= BASE_RADIUS / 4.0;
						ny *= BASE_RADIUS / 4.0;
						out
								.drawLine((orig.x + nx), (orig.y + ny),
										(dest.x + nx), (dest.y + ny),
										conf._bpThickness);
						out
								.drawLine((orig.x - nx), (orig.y - ny),
										(dest.x - nx), (dest.y - ny),
										conf._bpThickness);
					}
				} else if (style.isWobbleUG()) {
					double cx = (dest.x + orig.x) / 2.0;
					double cy = (dest.y + orig.y) / 2.0;
					out.drawLine(orig.x, orig.y, dest.x, dest.y,
							conf._bpThickness);
					drawSymbol(out, cx, cy, nx, ny, radiusCircle,
							style.isCIS(), style.getEdgePartner5(), thickness);
				} else {
					out.drawLine(orig.x, orig.y, dest.x, dest.y,
							conf._bpThickness);
				}
			} else {
				ModeleStyleBP.Edge p1 = style.getEdgePartner5();
				ModeleStyleBP.Edge p2 = style.getEdgePartner3();
				double cx = (dest.x + orig.x) / 2.0;
				double cy = (dest.y + orig.y) / 2.0;
				out.drawLine(orig.x, orig.y, dest.x, dest.y, conf._bpThickness);
				if (p1 == p2) {
					drawSymbol(out, cx, cy, nx, ny, radiusCircle,
							style.isCIS(), p1, thickness);
				} else {
					double vdx = (dest.x - orig.x);
					double vdy = (dest.y - orig.y);
					vdx /= 6.0;
					vdy /= 6.0;
					drawSymbol(out, cx + vdx, cy + vdy, nx, ny, radiusCircle,
							style.isCIS(), p2, thickness);
					drawSymbol(out, cx - vdx, cy - vdy, nx, ny, radiusCircle,
							style.isCIS(), p1, thickness);
				}
			}
		} else if (conf._mainBPStyle == VARNAConfig.BP_STYLE.BP_STYLE_RNAVIZ) {
			double xcenter = (orig.x + dest.x) / 2.0;
			double ycenter = (orig.y + dest.y) / 2.0;
			out.fillCircle(xcenter, ycenter, 3.0 * conf._bpThickness,
					conf._bpThickness, out.getCurrentColor());
		} else if (conf._mainBPStyle == VARNAConfig.BP_STYLE.BP_STYLE_SIMPLE) {
			out.drawLine(orig.x, orig.y, dest.x, dest.y, conf._bpThickness);
		}
	}

	
	private void drawColorMap(VARNAConfig _conf,SecStrDrawingProducer out)
	{
		double v1 = _conf._cm.getMinValue();
		double v2 = _conf._cm.getMaxValue();
		int x,y;
		double xSpaceAvail = 0;
		double ySpaceAvail = 0;
		double thickness = 1.0;
		/*ySpaceAvail = Math.min((getHeight()-rnabbox.height*scaleFactor-getTitleHeight())/2.0,scaleFactor*(_conf._colorMapHeight+VARNAConfig.DEFAULT_COLOR_MAP_FONT_SIZE));
		if ((int)ySpaceAvail==0)
		{
			xSpaceAvail = Math.min((getWidth()-rnabbox.width*scaleFactor)/2,scaleFactor*(_conf._colorMapWidth)+VARNAConfig.DEFAULT_COLOR_MAP_STRIPE_WIDTH);			
		}*/
		Rectangle2D.Double currentBBox = out.getBoundingBox(); 
		
		double xBase =  (currentBBox.getMaxX() -_conf._colorMapWidth-_conf._colorMapXOffset);
		//double yBase =  (minY - _conf._colorMapHeight + _conf._colorMapYOffset);
		double yBase =  (currentBBox.getMinY()-_conf._colorMapHeight-VARNAConfig.DEFAULT_COLOR_MAP_FONT_SIZE);
		
		for (int i=0;i<_conf._colorMapWidth;i++)
		{
			double ratio = (((double)i)/((double)_conf._colorMapWidth-1));
			double val = v1+(v2-v1)*ratio;
			Color c = _conf._cm.getColorForValue(val);
			x = (int) (xBase + i);
			y = (int) yBase;
			out.fillRectangle(x, y, VARNAConfig.DEFAULT_COLOR_MAP_STRIPE_WIDTH, _conf._colorMapHeight,c);	
		}
		out.setColor(VARNAConfig.DEFAULT_COLOR_MAP_OUTLINE);
		out.drawRectangle(xBase,yBase, (double)_conf._colorMapWidth+VARNAConfig.DEFAULT_COLOR_MAP_STRIPE_WIDTH-1, _conf._colorMapHeight,thickness);
		
		out.setColor(VARNAConfig.DEFAULT_COLOR_MAP_FONT_COLOR);
		out.setFont(out.getCurrentFont(),VARNAConfig.DEFAULT_COLOR_MAP_FONT_SIZE/1.5);
		out.drawText(xBase,
				yBase+_conf._colorMapHeight+VARNAConfig.DEFAULT_COLOR_MAP_FONT_SIZE/1.7,
				""+_conf._cm.getMinValue());
		out.drawText(xBase+VARNAConfig.DEFAULT_COLOR_MAP_STRIPE_WIDTH+_conf._colorMapWidth,
				yBase+_conf._colorMapHeight+VARNAConfig.DEFAULT_COLOR_MAP_FONT_SIZE/1.7,
				""+_conf._cm.getMaxValue());
			out.drawText(xBase+(VARNAConfig.DEFAULT_COLOR_MAP_STRIPE_WIDTH+_conf._colorMapWidth)/2.0,
				yBase-(VARNAConfig.DEFAULT_COLOR_MAP_FONT_SIZE/1.7),
				_conf._colorMapCaption);
		
	}

	
	private void renderRegionHighlights(SecStrDrawingProducer out, Point2D.Double[] realCoords,Point2D.Double[] realCenters)
	{
		for (HighlightRegionAnnotation r:_listeRegionHighlights)
		{
			GeneralPath s = r.getShape(realCoords,realCenters,1.0);
			out.setColor(r.getFillColor());
			out.fillPolygon(s, r.getFillColor());
			out.setColor(r.getOutlineColor());
			out.drawPolygon(s, 1l);
		}

	}
	
	private void saveRNA(String path, VARNAConfig conf, double scale,
			SecStrDrawingProducer out) throws ExceptionWritingForbidden {
		out.setScale(scale);
		// Computing bounding boxes
		double EPSMargin = 40;
		double minX = Double.MAX_VALUE;
		double maxX = Double.MIN_VALUE;
		double minY = Double.MAX_VALUE;
		double maxY = Double.MIN_VALUE;

		double x0, y0, x1, y1, xc, yc, xp, yp, dx, dy, norm;

		for (int i = 0; i < _listeBases.size(); i++) {
			minX = Math.min(minX, (_listeBases.get(i).getCoords().getX()
					- BASE_RADIUS - EPSMargin));
			minY = Math.min(minY, -(_listeBases.get(i).getCoords().getY()
					- BASE_RADIUS - EPSMargin));
			maxX = Math.max(maxX, (_listeBases.get(i).getCoords().getX()
					+ BASE_RADIUS + EPSMargin));
			maxY = Math.max(maxY, -(_listeBases.get(i).getCoords().getY()
					+ BASE_RADIUS + EPSMargin));
		}
		
		// Rescaling everything
		Point2D.Double[] coords = new Point2D.Double[_listeBases.size()]; 
		Point2D.Double[] centers = new Point2D.Double[_listeBases.size()]; 
		for (int i = 0; i < _listeBases.size(); i++) {
			xp = (_listeBases.get(i).getCoords().getX() - minX);
			yp = -(_listeBases.get(i).getCoords().getY() - minY);
			coords[i] = new Point2D.Double(xp,yp);

			Point2D.Double centerBck = getCenter(i);
			if (get_drawMode() == RNA.DRAW_MODE_NAVIEW
					|| get_drawMode() == RNA.DRAW_MODE_RADIATE) 
			{
				if ((_listeBases.get(i).getElementStructure() != -1)
						&& i < _listeBases.size() - 1
						&& i > 1)
				{
				  ModeleBase b1 = get_listeBases().get(i - 1);
				  ModeleBase b2 = get_listeBases().get(i + 1);
				  int j1 = b1.getElementStructure();
				  int j2 = b2.getElementStructure();
				  if ((j1==-1)^(j2==-1))
				  {
				  // alors la position du nombre associé doit etre décalé
					Point2D.Double a1 = b1.getCoords();
					Point2D.Double a2 = b2.getCoords();
					Point2D.Double c1 = b1.getCenter();
					Point2D.Double c2 = b2.getCenter();
					
					centerBck.x = _listeBases.get(i).getCoords().x + (c1.x-a1.x)/c1.distance(a1)+(c2.x-a2.x)/c2.distance(a2);
					centerBck.y = _listeBases.get(i).getCoords().y + (c1.y-a1.y)/c1.distance(a1)+(c2.y-a2.y)/c2.distance(a2);
				  }
				}
			}
			xc = (centerBck.getX() - minX);
			yc = -(centerBck.getY() - minY);
			centers[i] = new Point2D.Double(xc,yc);
		}


		// Drawing background
		if (conf._drawBackground)
		  out.setBackgroundColor(conf._backgroundColor);
		
		// Drawing region highlights
		renderRegionHighlights(out,coords,centers);
		
		// Drawing backbone
		out.setColor(conf._backboneColor);
		for (int i = 1; i < _listeBases.size(); i++) {
			x0 = coords[i-1].x;
			y0 = coords[i-1].y;
			x1 = coords[i].x;
			y1 = coords[i].y;
			Point2D.Double vn = new Point2D.Double();
			double dist = coords[i-1].distance(coords[i]);
			boolean discontinuous= (getBaseNumber(i)-getBaseNumber(i-1)!=1);
			if ((dist>0)&!(discontinuous)){
					vn.x = (x1-x0)/dist;
					vn.y = (y1-y0)/dist;
					out.drawLine((x0+BASE_RADIUS*vn.x),
						(y0+BASE_RADIUS*vn.y),
						(x1-BASE_RADIUS*vn.x),
						(y1-BASE_RADIUS*vn.y),1.0);
			}
		}

		// Drawing bonds
		for (int i = 0; i < _listeBases.size(); i++) {
			if (_listeBases.get(i).getElementStructure() > i) {
				ModeleStyleBP style = _listeBases.get(i).getStyleBP();
				if (style.isCanonical() || conf._drawnNonCanonicalBP) {
					Color bpcol = getBasePairColor(style, conf);
					out.setColor(bpcol);

					int j = _listeBases.get(i).getElementStructure();
					x0 = coords[i].x;
					y0 = coords[i].y;
					x1 = coords[j].x;
					y1 = coords[j].y;
					dx = x1 - x0;
					dy = y1 - y0;
					norm = Math.sqrt(dx * dx + dy * dy);
					dx /= norm;
					dy /= norm;

					if (_drawMode == DRAW_MODE_CIRCULAR
							|| _drawMode == DRAW_MODE_RADIATE
							|| _drawMode == DRAW_MODE_NAVIEW) {
						drawBasePair(out, new Point2D.Double(x0, y0),
								new Point2D.Double(x1, y1), style, conf);
					} else if (_drawMode == DRAW_MODE_LINEAR) {
						double coef;
						double distance;
						if (j - i == 1)
							coef = _bpHeightIncrement * 2;
						else
							coef = _bpHeightIncrement * 1;
						distance = (int) Math.round(x1 - x0);
						out.drawArc(new Point2D.Double(x0, y0), distance,
								distance * coef, 0, 180);
					}
				}
			}
		}

		// Drawing additional bonds
		if (conf._drawnNonPlanarBP) {
			for (int i = 0; i < _structureAux.size(); i++) {
				ModeleStyleBP bp = _structureAux.get(i);
				out.setColor(getBasePairColor(bp, conf));
				
				int a = bp.getPartner5().getIndex();
				int b = bp.getPartner3().getIndex();

				if (bp.isCanonical() || conf._drawnNonCanonicalBP) {
					x0 = coords[a].x;
					y0 = coords[a].y;
					x1 = coords[b].x;
					y1 = coords[b].y;
					dx = x1 - x0;
					dy = y1 - y0;
					norm = Math.sqrt(dx * dx + dy * dy);
					dx /= norm;
					dy /= norm;
					if ((_drawMode == DRAW_MODE_CIRCULAR)
							|| (_drawMode == DRAW_MODE_RADIATE)
							|| _drawMode == DRAW_MODE_NAVIEW) {
						drawBasePair(out, new Point2D.Double(x0, y0),
								new Point2D.Double(x1, y1), bp, conf);
					} else if (_drawMode == DRAW_MODE_LINEAR) {
						double coef;
						double distance;
						if (b - a == 1)
							coef = _bpHeightIncrement * 2;
						else
							coef = _bpHeightIncrement * 1;
						distance = (int) Math.round(x1 - x0);
						out.drawArc(new Point2D.Double(x0, y0), distance,
								distance * coef, 0, 180);
					}
				}
			}
		}

		// Drawing Bases
		double baseFontSize = (1.5 * BASE_RADIUS);
		out.setFont(PSExport.FONT_HELVETICA_BOLD, baseFontSize);
		if (_comparisonMode) {
			for (int i = 0; i < _listeBases.size(); i++) {
				x0 = coords[i].x;
				y0 = coords[i].y;
				out.fillCircle(x0, y0, BASE_RADIUS, 1l, getBaseInnerColor(i,
						conf));
				out.setColor(getBaseOuterColor(i, conf));
				out.drawCircle(x0, y0, BASE_RADIUS, 1l);
				out.setColor(getBaseNameColor(i, conf));
				out.drawText(x0, y0, ((ModeleBasesComparison) _listeBases
						.get(i)).getBases());
			}
		} else {
			for (int i = 0; i < _listeBases.size(); i++) {
				ModeleBase mb = _listeBases.get(i);
				x0 = coords[i].x;
				y0 = coords[i].y;
				if (conf._fillBase)
				{
				  out.fillCircle(x0, y0, BASE_RADIUS, 1l, getBaseInnerColor(i,
						conf));
				}
				if (conf._drawOutlineBase)
				{
				out.setColor(getBaseOuterColor(i, conf));
				out.drawCircle(x0, y0, BASE_RADIUS, 1l);
				}
				out.setColor(getBaseNameColor(i, conf));
				out.drawText(x0, y0, _listeBases.get(i).getContent());
			}
		}

		// Drawing base numbers
		double numFontSize = (double) (1.5 * BASE_RADIUS);
		out.setFont(PSExport.FONT_HELVETICA_BOLD, numFontSize);

		for (int i = 0; i < _listeBases.size(); i++) {
			int basenum = _listeBases.get(i).getBaseNumber();
			if (basenum == -1) {
				basenum = i + 1;
			}
			if (this.isNumberDrawn(_listeBases.get(i),conf._numPeriod)) {
				out.setColor(_listeBases.get(i).getStyleBase().get_base_number_color());
				x0 = coords[i].x;
				y0 = coords[i].y;
				x1 = centers[i].x;
				y1 = centers[i].y;
				dx = x1 - x0;
				dy = y1 - y0;
				norm = Math.sqrt(dx * dx + dy * dy);
				dx /= norm;
				dy /= norm;
				out.drawLine((x0 - 1.5 * BASE_RADIUS * dx), (y0 - 1.5
						* BASE_RADIUS * dy), (x0 - 2.5 * BASE_RADIUS * dx),
						(y0 - 2.5 * BASE_RADIUS * dy), 1);
				out.drawText((x0 - (conf._distNumbers+1.0) * BASE_RADIUS * dx), (y0 - (conf._distNumbers+1.0) * BASE_RADIUS
						* dy), "" + (basenum));
			}
		}
		renderAnnotations(out, minX, minY,conf);
		
		// Draw color map 
		if (conf._drawColorMap)
		{ drawColorMap(conf, out); }

		
		// Drawing Title
		Rectangle2D.Double currentBBox = out.getBoundingBox(); 
		double titleFontSize = (2.0*conf._titleFont.getSize());
		out.setColor(conf._titleColor);
		out.setFont(PSExport.FONT_HELVETICA_BOLD,titleFontSize);
		double yTitle = currentBBox.y-titleFontSize/2.0;
		if (!getName().equals(""))
		{
			out.drawText((maxX-minX)/2.0,yTitle, getName());			
		}
		else if (!conf._title.equals(""))
		{
			out.drawText((maxX-minX)/2.0,yTitle, getName());			
		}


		FileWriter fout;
		try {
			fout = new FileWriter(path);
			fout.write(out.export());
			fout.close();
		} catch (IOException e) {
			throw new ExceptionWritingForbidden(e.getMessage());
		}
	}
	
	Point2D.Double buildCaptionPosition(ModeleBase mb, double heightEstimate, VARNAConfig conf)
	{
		double radius = 2.0;
		if (isNumberDrawn(mb,conf._numPeriod))
		{ radius += (conf._distNumbers+1.0); }
		Point2D.Double center = mb.getCenter();
		Point2D.Double p = mb.getCoords();
		double realDistance = BASE_RADIUS*radius + heightEstimate;
		return new Point2D.Double(center.getX()
		+ (p.getX() - center.getX())
		* ((p.distance(center) + realDistance) / p
				.distance(center)), center.getY()
		+ (p.getY() - center.getY())
		* ((p.distance(center) + realDistance) / p
				.distance(center)) );
	}

    public double getBPHeightIncrement()
    {
    	return this._bpHeightIncrement;
    }
	
    public void setBPHeightIncrement(double d)
    {
    	_bpHeightIncrement = d;
    }
	
	public static double CHEM_PROB_ARROW_THICKNESS = 2.0;
	
	private void drawChemProbAnnotation(SecStrDrawingProducer out, ChemProbAnnotation cpa, Point2D.Double anchor, double minX, double minY)
	{
		out.setColor(cpa.getColor());
		Point2D.Double v = cpa.getDirVector();
		Point2D.Double vn = cpa.getNormalVector();
		Point2D.Double base = new Point2D.Double((anchor.x+CHEM_PROB_DIST*v.x),(anchor.y+CHEM_PROB_DIST*v.y));
		Point2D.Double edge = new Point2D.Double((base.x+CHEM_PROB_BASE_LENGTH*cpa.getIntensity()*v.x),(base.y+CHEM_PROB_BASE_LENGTH*cpa.getIntensity()*v.y));
		double thickness = CHEM_PROB_ARROW_THICKNESS*cpa.getIntensity();
		switch (cpa.getType())
		{
		  case ARROW_TYPE:
		  {
			  Point2D.Double arrowTip1 = new Point2D.Double((base.x+cpa.getIntensity()*(CHEM_PROB_ARROW_WIDTH*vn.x+CHEM_PROB_ARROW_HEIGHT*v.x)),
					  (base.y+cpa.getIntensity()*(CHEM_PROB_ARROW_WIDTH*vn.y+CHEM_PROB_ARROW_HEIGHT*v.y)));
			  Point2D.Double arrowTip2 = new Point2D.Double((base.x+cpa.getIntensity()*(-CHEM_PROB_ARROW_WIDTH*vn.x+CHEM_PROB_ARROW_HEIGHT*v.x)),
					  (base.y+cpa.getIntensity()*(-CHEM_PROB_ARROW_WIDTH*vn.y+CHEM_PROB_ARROW_HEIGHT*v.y)));
			  out.drawLine(base.x-minX,minY-base.y,edge.x-minX,minY-edge.y,thickness);
			  out.drawLine(base.x-minX,minY-base.y,arrowTip1.x-minX,minY-arrowTip1.y,thickness);
			  out.drawLine(base.x-minX,minY-base.y,arrowTip2.x-minX,minY-arrowTip2.y,thickness);
		  }
		  break;
		  case PIN_TYPE:
		  {
			  Point2D.Double side1 = new Point2D.Double((edge.x-cpa.getIntensity()*(CHEM_PROB_PIN_SEMIDIAG*v.x)),
					  (edge.y-cpa.getIntensity()*(CHEM_PROB_PIN_SEMIDIAG*v.y)));
			  Point2D.Double side2 = new Point2D.Double((edge.x-cpa.getIntensity()*(CHEM_PROB_PIN_SEMIDIAG*vn.x)),
					  (edge.y-cpa.getIntensity()*(CHEM_PROB_PIN_SEMIDIAG*vn.y)));
			  Point2D.Double side3 = new Point2D.Double((edge.x+cpa.getIntensity()*(CHEM_PROB_PIN_SEMIDIAG*v.x)),
					  (edge.y+cpa.getIntensity()*(CHEM_PROB_PIN_SEMIDIAG*v.y)));
			  Point2D.Double side4 = new Point2D.Double((edge.x+cpa.getIntensity()*(CHEM_PROB_PIN_SEMIDIAG*vn.x)),
					  (edge.y+cpa.getIntensity()*(CHEM_PROB_PIN_SEMIDIAG*vn.y)));
				GeneralPath p2 = new GeneralPath();
				p2.moveTo((float)(side1.x-minX),(float)(minY-side1.y));
				p2.lineTo((float)(side2.x-minX),(float)(minY-side2.y));
				p2.lineTo((float)(side3.x-minX),(float)(minY-side3.y));
				p2.lineTo((float)(side4.x-minX),(float)(minY-side4.y));
				p2.closePath();
				out.fillPolygon(p2, cpa.getColor());
				out.drawLine(base.x-minX,minY-base.y,edge.x-minX,minY-edge.y,thickness);
		  }
		  break;
		  case TRIANGLE_TYPE:
		  {
			  Point2D.Double arrowTip1 = new Point2D.Double((edge.x+cpa.getIntensity()*(CHEM_PROB_TRIANGLE_WIDTH*vn.x)),
					  (edge.y+cpa.getIntensity()*(CHEM_PROB_TRIANGLE_WIDTH*vn.y)));
			  Point2D.Double arrowTip2 = new Point2D.Double((edge.x+cpa.getIntensity()*(-CHEM_PROB_TRIANGLE_WIDTH*vn.x)),
					  (edge.y+cpa.getIntensity()*(-CHEM_PROB_TRIANGLE_WIDTH*vn.y)));
				GeneralPath p2 = new GeneralPath();
				p2.moveTo((float)(base.x-minX),(float)(minY-base.y));
				p2.lineTo((float)(arrowTip1.x-minX),(float)(minY-arrowTip1.y));
				p2.lineTo((float)(arrowTip2.x-minX),(float)(minY-arrowTip2.y));
				p2.closePath();
				out.fillPolygon(p2, cpa.getColor());
		  }
		  break;
		  case DOT_TYPE:
		  {
			  Double radius = CHEM_PROB_DOT_RADIUS*cpa.getIntensity();
			  Point2D.Double center = new Point2D.Double((base.x+radius*v.x)-minX,minY-(base.y+radius*v.y));
			  out.fillCircle(center.x, center.y, radius, thickness, cpa.getColor());
		  }
		  break;
		}
	}

	
	private void renderAnnotations(SecStrDrawingProducer out, double minX, double minY, VARNAConfig conf) 
	{
		for (TextAnnotation textAnnotation : getAnnotations()) {
			out.setColor(textAnnotation.getColor());
			out.setFont(PSExport.FONT_HELVETICA_BOLD,2.0*textAnnotation.getFont().getSize());
			Point2D.Double position = textAnnotation.getCenterPosition();
			if (textAnnotation.getType()==TextAnnotation.BASE)
			{
				ModeleBase mb = (ModeleBase) textAnnotation.getAncrage();
				double fontHeight = Math.ceil(textAnnotation.getFont().getSize());
				position = buildCaptionPosition(mb,fontHeight,conf);

			}
			out.drawText(position.x-minX,-(position.y-minY), textAnnotation.getTexte());
		}
		for (ChemProbAnnotation cpa : getChemProbAnnotations()) {
			Point2D.Double anchor = cpa.getAnchorPosition();
			drawChemProbAnnotation(out,cpa,anchor,minX,minY);
		}
	}


	public boolean isNumberDrawn(ModeleBase mb, int numPeriod)
	{
		return (mb.getIndex() == 0) || ((mb.getBaseNumber()) % numPeriod == 0)
		|| (mb.getIndex() == get_listeBases().size() - 1);
	}
	
	
	public void saveRNAEPS(String path, VARNAConfig conf)
			throws ExceptionWritingForbidden {
		PSExport out = new PSExport();
		saveRNA(path, conf, 0.4, out);
	}

	public void saveRNAXFIG(String path, VARNAConfig conf)
			throws ExceptionWritingForbidden {
		XFIGExport out = new XFIGExport();
		saveRNA(path, conf, 20, out);
	}

	public void saveRNASVG(String path, VARNAConfig conf)
			throws ExceptionWritingForbidden {
		SVGExport out = new SVGExport();
		saveRNA(path, conf, 0.5, out);
	}

	public Rectangle2D.Double getBBox() {
		Rectangle2D.Double result = new Rectangle2D.Double(10, 10, 10, 10);
		double minx, maxx, miny, maxy;
		minx = Double.MAX_VALUE;
		miny = Double.MAX_VALUE;
		maxx = -Double.MAX_VALUE;
		maxy = -Double.MAX_VALUE;
		for (int i = 0; i < _listeBases.size(); i++) {
			minx = Math.min(_listeBases.get(i).getCoords().getX()
					- BASE_RADIUS, minx);
			miny = Math.min(_listeBases.get(i).getCoords().getY()
					- BASE_RADIUS, miny);
			maxx = Math.max(_listeBases.get(i).getCoords().getX()
					+ BASE_RADIUS, maxx);
			maxy = Math.max(_listeBases.get(i).getCoords().getY()
					+ BASE_RADIUS, maxy);
		}
		result.x = minx;
		result.y = miny;
		result.width = Math.max(maxx - minx,1);
		result.height = Math.max(maxy - miny,1);
		if (_drawMode==RNA.DRAW_MODE_LINEAR)
		{
			double realHeight = _bpHeightIncrement*result.width/2.0;
			result.height += realHeight;
			result.y -= realHeight; 
		}
		return result;
	}

	public void setCoord(int index, Point2D.Double p) {
		setCoord(index, p.x, p.y);
	}

	public void setCoord(int index, double x, double y) {
		if (index < _listeBases.size()) {
			_listeBases.get(index).setCoords(new Point2D.Double(x, y));
		}
	}

	public Point2D.Double getCoords(int i) {
		if (i < _listeBases.size() && i >= 0) {
			return _listeBases.get(i).getCoords();
		}
		return new Point2D.Double();
	}

	public String getBaseContent(int i) {
		if ((i>=0)&&(i < _listeBases.size())) 
		{ return _listeBases.get(i).getContent(); }
		return "";
	}

	public int getBaseNumber(int i) {
		if ((i>=0)&&(i < _listeBases.size())) 
		{ return _listeBases.get(i).getBaseNumber(); }
		return -1;
	}
	
	public Point2D.Double getCenter(int i) {
		if (i < _listeBases.size()) {
			return _listeBases.get(i).getCenter();
		}

		return new Point2D.Double();
	}

	public void setCenter(int i, Point2D.Double p) {
		if (i < _listeBases.size()) {
			_listeBases.get(i).setCenter(p);
		}
	}

	public void drawRNACircle() {
		_drawn = true;
		_drawMode = DRAW_MODE_CIRCULAR;
		int radius = (int) ((3 * (_listeBases.size() + 1) * BASE_RADIUS) / (2 * Math.PI));
		double angle;
		for (int i = 0; i < _listeBases.size(); i++) {
			angle = -((((double) -(i + 1)) * 2.0 * Math.PI)
					/ ((double) (_listeBases.size() + 1)) - Math.PI / 2.0);
			_listeBases.get(i).setCoords(
					new Point2D.Double(
							(radius * Math.cos(angle) * _spaceBetweenBases),
							(radius * Math.sin(angle) * _spaceBetweenBases)));
			_listeBases.get(i).setCenter(new Point2D.Double(0, 0));
		}
	}

	public void drawRNAVARNAView() {
		_drawn = true;
		_drawMode = DRAW_MODE_VARNA_VIEW;
		VARNASecDraw vs = new VARNASecDraw();
		vs.drawRNA(1, this);
	}
	
	
	public void drawRNALine() {
		_drawn = true;
		_drawMode = DRAW_MODE_LINEAR;
		for (int i = 0; i < get_listeBases().size(); i++) {
			get_listeBases().get(i).setCoords(
					new Point2D.Double(i * _spaceBetweenBases * 20, 0));
			get_listeBases().get(i).setCenter(
					new Point2D.Double(i * _spaceBetweenBases * 20, -10));
		}
	}
	
	
	/**
	 * Argument helixEndPoint is an IN argument (will be read),
	 * and must contain an helix edge endpoint.
	 * 
	 * Arguments helixEndPointCoords, helixVector and j are OUT arguments
	 * (must be existing objects, content will be overwritten).
	 * 
	 * The helixEndPointCoords argument will contain the coordinates of
	 * this helix endpoint.
	 * 
	 * The helixVector argument will contain the vector
	 * from the helix startPosition to endPosition or the opposite
	 * depending on there the endpoint is (the endpoint will be on the
	 * destination side of the vector).
	 * 
	 * The j vector will contain an unitary vector that is colinear
	 * to the last/first base pair connection on the side of this endpoint.
	 * The vector will be oriented to the side of the given endpoint.
	 */
	private void whereIsThisHelixEndpoint(
			EdgeEndPoint helixEndPoint,
			Point2D.Double helixEndPointCoords,
			Point2D.Double helixVector,
			Point2D.Double j) {
		RNATemplateHelix helix = (RNATemplateHelix) helixEndPoint.getElement();
		Point2D.Double startpos = helix.getStartPosition();
		Point2D.Double endpos = helix.getEndPosition();
		Point2D.Double thisSideHelixPos = new Point2D.Double();
		switch (helixEndPoint.getPosition()) {
		case IN1:
		case OUT2:
			helixVector.x = startpos.x - endpos.x;
			helixVector.y = startpos.y - endpos.y;
			thisSideHelixPos = startpos;
			break;
		case IN2:
		case OUT1:
			helixVector.x = endpos.x - startpos.x;
			helixVector.y = endpos.y - startpos.y;
			thisSideHelixPos = endpos;
			break;
		}
		double helixVectorLength = Math.hypot(helixVector.x, helixVector.y);
		// i is the vector which is colinear to helixVector and such that ||i|| = 1
		Point2D.Double i = new Point2D.Double();
		i.x = helixVector.x / helixVectorLength;
		i.y = helixVector.y / helixVectorLength;
		// Find j such that it is orthogonal to i, ||j|| = 1
		// and j goes to the side where the sequence will be connected
		switch (helixEndPoint.getPosition()) {
		case IN1:
		case IN2:
			// rotation of +pi/2
			j.x = - i.y;
			j.y =   i.x;
			break;
		case OUT1:
		case OUT2:
			// rotation of -pi/2
			j.x =   i.y;
			j.y = - i.x;
			break;
		}
		if (helix.isFlipped()) {
			j.x = - j.x;
			j.y = - j.y;
		}
		// Compute the position of this helix endpoint
		helixEndPointCoords.x = thisSideHelixPos.x + BASE_PAIR_DISTANCE / 2 * j.x;
		helixEndPointCoords.y = thisSideHelixPos.y + BASE_PAIR_DISTANCE / 2 * j.y;
	}
	
	/**
	 * A cubic Bezier curve can be defined by 4 points,
	 * see http://en.wikipedia.org/wiki/Bezier_curve#Cubic_B.C3.A9zier_curves
	 * For each of the curve end points, there is the last/first point of the
	 * curve and a point which gives the direction and length of the tangent
	 * vector on that side. This two points are respectively curveEndPoint
	 * and curveVectorOtherPoint.
	 * Argument helixVector is the vector formed by the helix,
	 * in the right direction for our sequence.
	 * All this Point2D.Double are "OUT arguments".
	 * They must be allocated and the values will be modified.
	 * Argument sequenceEndPointIsIn is true if we consider the "in" endpoint
	 * of the sequence, and false if we consider the "out" endpoint.
	 */
	private void computeBezierPoints(EdgeEndPoint helixEndPoint,
			Point2D.Double curveEndPoint,
			Point2D.Double curveVectorOtherPoint,
			Point2D.Double helixVector)
			throws RNATemplateDrawingAlgorithmException {
		
		RNATemplateUnpairedSequence sequence = (RNATemplateUnpairedSequence) helixEndPoint.getOtherElement();
		EdgeEndPointPosition endPointPositionOnHelix = helixEndPoint.getPosition();
		boolean sequenceEndPointIsIn;
		switch (endPointPositionOnHelix) {
		case IN1:
		case IN2:
			sequenceEndPointIsIn = false;
			break;
		default:
			sequenceEndPointIsIn = true;
		}
		
		EdgeEndPoint endPointOnHelix =
			sequenceEndPointIsIn ?
					sequence.getIn().getOtherEndPoint() :
					sequence.getOut().getOtherEndPoint();
		if (endPointOnHelix == null) {
			throw (new RNATemplateDrawingAlgorithmException("Sequence is not connected to an helix."));
		}
		RNATemplateHelix helix = (RNATemplateHelix) endPointOnHelix.getElement();
		
		double l =
			sequenceEndPointIsIn ?
				sequence.getInTangentVectorLength() :
				sequence.getOutTangentVectorLength();
		Point2D.Double startpos = helix.getStartPosition();
		Point2D.Double endpos = helix.getEndPosition();
		Point2D.Double thisSideHelixPos = new Point2D.Double();
		if ((helix.getOut1() == endPointOnHelix && sequenceEndPointIsIn)
				|| (helix.getIn2() == endPointOnHelix && !sequenceEndPointIsIn)) {
			helixVector.x = endpos.x - startpos.x;
			helixVector.y = endpos.y - startpos.y;
			thisSideHelixPos = endpos;
		} else if ((helix.getOut2() == endPointOnHelix && sequenceEndPointIsIn)
				|| (helix.getIn1() == endPointOnHelix && !sequenceEndPointIsIn)) {
			helixVector.x = startpos.x - endpos.x;
			helixVector.y = startpos.y - endpos.y;
			thisSideHelixPos = startpos;
		} else {
			throw (new RNATemplateDrawingAlgorithmException("Connection problem between helix and sequence."));
		}
		double helixVectorLength = Math.hypot(helixVector.x, helixVector.y);
		// i is the vector which is colinear to helixVector and such that ||i|| = 1
		Point2D.Double i = new Point2D.Double();
		i.x = helixVector.x / helixVectorLength;
		i.y = helixVector.y / helixVectorLength;
		// Find j such that it is orthogonal to i, ||j|| = 1
		// and j goes to the side where the sequence will be connected
		Point2D.Double j = new Point2D.Double();
		if (!sequenceEndPointIsIn) {
			j.x = - i.y;
			j.y =   i.x;
		} else {
			j.x =   i.y;
			j.y = - i.x;
		}
		if (helix.isFlipped()) {
			j.x = - j.x;
			j.y = - j.y;
		}
		// Compute the endpoint of the curve
		curveEndPoint.x = thisSideHelixPos.x + BASE_PAIR_DISTANCE / 2 * j.x;
		curveEndPoint.y = thisSideHelixPos.y + BASE_PAIR_DISTANCE / 2 * j.y;
		
		// Compute the absolute angle our line makes to the helix
		double theta =
			sequenceEndPointIsIn ?
				sequence.getInTangentVectorAngle() :
				sequence.getOutTangentVectorAngle();
		
		// Compute v, the tangent vector of the Bezier curve
		Point2D.Double v = new Point2D.Double();
		v.x = l * Math.cos(theta);
		v.y = l * Math.sin(theta);
		curveVectorOtherPoint.x = curveEndPoint.x + v.x;
		curveVectorOtherPoint.y = curveEndPoint.y + v.y;
	}
	

	/**
	 * Compute the angle made by a vector.
	 */
	private static double angleFromVector(Point2D.Double v) {
		double l = Math.hypot(v.x, v.y);
		if (v.y > 0) {
			return Math.acos(v.x / l);
		} else if (v.y < 0) {
			return - Math.acos(v.x / l);
		} else {
			return v.x > 0 ? 0 : Math.PI;
		}
	}
	
	private static double angleFromVector(double x, double y) {
		return angleFromVector(new Point2D.Double(x, y));
	}
	
	/**
	 * Draw the given helix (given as a *sorted* array of indexes)
	 * like defined in the given template helix.
	 * The bases positions are not changed in fact, instead the coords and
	 * centers arrays are modified.
	 */
	private void drawHelixLikeTemplateHelix(int[] basesInHelixArray,
			RNATemplateHelix helix,
			Point2D.Double[] coords,
			Point2D.Double[] centers) {
		int n = basesInHelixArray.length / 2;
		Point2D.Double startpos, endpos;
		if (helix.getIn1Is() == In1Is.IN1_IS_5PRIME) {
			startpos = helix.getStartPosition();
			endpos = helix.getEndPosition();
		} else {
			endpos = helix.getStartPosition();
			startpos = helix.getEndPosition();
		}
		// (i_x,i_y) is the vector between to consecutive bases
		// of the same side of an helix
		Point2D.Double i = new Point2D.Double();
		i.x = (endpos.x - startpos.x) / (n-1);
		i.y = (endpos.y - startpos.y) / (n-1);
		Point2D.Double j = new Point2D.Double();
		j.x = - i.y;
		j.y =   i.x;
		if (helix.isFlipped()) {
			j.x = - j.x;
			j.y = - j.y;
		}
		double j_original_norm = Math.hypot(j.x, j.y);
		// change (j_x,j_y) so that its norm is 1
		j.x = j.x / j_original_norm;
		j.y = j.y / j_original_norm;
		Point2D.Double o = new Point2D.Double();
		o.x = startpos.x - j.x * BASE_PAIR_DISTANCE / 2;
		o.y = startpos.y - j.y * BASE_PAIR_DISTANCE / 2;
		
		for (int k=0; k<n; k++) {
			int b1 = basesInHelixArray[k];
			int b2 = basesInHelixArray[2*n-k-1];
			Point2D.Double p1 = new Point2D.Double();
			p1.x = o.x + k*i.x;
			p1.y = o.y + k*i.y;
			coords[b1] = p1;
			Point2D.Double p2 = new Point2D.Double();
			p2.x = p1.x + j.x * BASE_PAIR_DISTANCE;
			p2.y = p1.y + j.y * BASE_PAIR_DISTANCE;
			coords[b2] = p2;
		}
	
		for (int k=0; k<2*n-1; k++) {
			if (k == n-1) continue;
			int b1 = basesInHelixArray[k];
			int b2 = basesInHelixArray[k+1];
			if (b1 + 1 != b2) {
				// There is a loop between these 2 bases
				Point2D.Double b1pos = coords[b1];
				Point2D.Double b2pos = coords[b2];
				// ||v|| = 1
				Point2D.Double v = new Point2D.Double();
				if (k >= n) {
					v.x = j.x;
					v.y = j.y;
				} else {
					v.x = - j.x;
					v.y = - j.y;
				}
				drawLoop(b1+1,
						 b2-1,
						 (b1pos.x + b2pos.x)/2 + v.x * LOOP_DISTANCE,
						 (b1pos.y + b2pos.y)/2 + v.y * LOOP_DISTANCE,
						 angleFromVector(v),
						 coords,
						 centers);
			}
		}
	}
	
	/**
	 * A Bezier curve can be defined by four points,
	 * see http://en.wikipedia.org/wiki/Bezier_curve#Cubic_B.C3.A9zier_curves
	 * Here we give this four points and an array of bases indexes
	 * (which must be indexes in this RNA sequence) which will be moved
	 * to be on the Bezier curve.
	 * The bases positions are not changed in fact, instead the coords and
	 * centers arrays are modified.
	 */
	private void drawOnBezierCurve(int[] basesInSequence,
			Point2D.Double P0,
			Point2D.Double P1,
			Point2D.Double P2,
			Point2D.Double P3,
			Point2D.Double[] coords,
			Point2D.Double[] centers) {
		// Draw the bases of the sequence along a Bezier curve
		int n = basesInSequence.length;
		// We choose to approximate the Bezier curve by 10*n straight lines.
		CubicBezierCurve bezier = new CubicBezierCurve(P0, P1, P2, P3, 10*n);
		double curveLength = bezier.getApproxCurveLength();
		double delta_t = curveLength / (n+1);
		double[] t = new double[n];
		for (int k=0; k<n; k++) {
			t[k] = (k+1) * delta_t;
		}
		Point2D.Double[] sequenceBasesCoords = bezier.uniformParam(t);
		for (int k=0; k<n; k++) {
			coords[basesInSequence[k]] = sequenceBasesCoords[k];
		}
	}
	
	/**
	 * Like drawOnBezierCurve(), but on a straight line.
	 */
	private void drawOnStraightLine(int[] basesInSequence,
			Point2D.Double P0,
			Point2D.Double P3,
			Point2D.Double[] coords,
			Point2D.Double[] centers) {
		// Draw the bases of the sequence along a Bezier curve
		int n = basesInSequence.length;
		Point2D.Double v = new Point2D.Double(P3.x - P0.x, P3.y - P0.y);
		for (int k=0; k<n; k++) {
			coords[basesInSequence[k]].x = P0.x + (k+1) / (double) (n+1) * v.x;
			coords[basesInSequence[k]].y = P0.y + (k+1) / (double) (n+1) * v.y;
		}
	}
	
	/**
	 * A Bezier curve can be defined by four points,
	 * see http://en.wikipedia.org/wiki/Bezier_curve#Cubic_B.C3.A9zier_curves
	 * This functions draws the RNA sequence between (including)
	 * firstBase and lastBase along the Bezier curve defined by (P0,P1,P2,P3).
	 * The sequence may contain helixes.
	 * The bases positions are not changed in fact, instead the coords and
	 * centers arrays are modified.
	 * If P1 and P2 are null, the bases are drawn on a straight line.
	 */
	private void drawAlongBezierCurve(
			int firstBase,
			int lastBase,
			Point2D.Double P0,
			Point2D.Double P1,
			Point2D.Double P2,
			Point2D.Double P3,
			Point2D.Double[] coords,
			Point2D.Double[] centers) {
		
		// First we find the bases which are directly on the Bezier curve
		ArrayList<Integer> alongBezierCurve = new ArrayList<Integer>();
		for (int depth=0, i=firstBase; i<=lastBase; i++) {
			int k = _listeBases.get(i).getElementStructure();
			if (k < 0) {
				if (depth == 0) {
					alongBezierCurve.add(i);
				}
			} else {
				if (i < k) {
					if (depth == 0) {
						alongBezierCurve.add(i);
						alongBezierCurve.add(k);
					}
					depth++;
				} else {
					depth--;
				}
			}
		}
		// number of bases along the Bezier curve
		int n = alongBezierCurve.size();
		int[] alongBezierCurveArray = RNATemplateAlign.intArrayFromList(alongBezierCurve);
		if (P1 != null && P2 != null) {
			drawOnBezierCurve(alongBezierCurveArray, P0, P1, P2, P3, coords, centers);
		} else {
			drawOnStraightLine(alongBezierCurveArray, P0, P3, coords, centers);
		}
		
		// Now use the radiate algorithm to draw the helixes
		for (int k=0; k<n-1; k++) {
			int b1 = alongBezierCurveArray[k];
			int b2 = alongBezierCurveArray[k+1];
			if (_listeBases.get(b1).getElementStructure() == b2) {
				Point2D.Double b1pos = coords[b1];
				Point2D.Double b2pos = coords[b2];
				double alpha = angleFromVector(b2pos.x - b1pos.x, b2pos.y - b1pos.y);
				drawLoop(b1,
						 b2,
						 (b1pos.x + b2pos.x)/2,
						 (b1pos.y + b2pos.y)/2,
						 alpha - Math.PI / 2,
						 coords,
						 centers);
			}
		}	
	}
	
	
	/**
	 * Starts at helixEndPoint and draws everything remaining until
	 * the sequence beginning or end, using the radiate algorithm.
	 * The helixEndPoint argument gives the base of the helix which is on
	 * the endpoint.
	 */
	private void drawFromHelixEndPoint(
			EdgeEndPoint helixEndPoint,
			int endPointLastBase,
			Point2D.Double[] coords,
			Point2D.Double[] centers) {
		Point2D.Double helixEndPointCoords = new Point2D.Double();
		Point2D.Double helixVector = new Point2D.Double();
		Point2D.Double fromEndPointVector = new Point2D.Double();
		RNATemplateHelix helix = (RNATemplateHelix) helixEndPoint.getElement();
		whereIsThisHelixEndpoint(helixEndPoint, helixEndPointCoords, helixVector, fromEndPointVector);
		
		// Find whether we are going in the 5'->3' direction or the 3'->5' direction
		int nextBase;
		switch (helixEndPoint.getPosition()) {
		case IN1:
		case IN2:
			nextBase = -1;
			break;
		default:
			nextBase = +1;
		}

		// This code was taken from the radiate algorithm
		double dirAngle = angleFromVector(fromEndPointVector);
		int i = endPointLastBase;
		double x = helixEndPointCoords.x;
		double y = helixEndPointCoords.y;
		double vx = -Math.sin(dirAngle);
		double vy = Math.cos(dirAngle);
		int baseCount = _listeBases.size();
		while(i >= 0 && i < baseCount)
		{
			coords[i].x = x; 
			coords[i].y = y;
			centers[i].x = x+BASE_PAIR_DISTANCE*vy;
			centers[i].y = y-BASE_PAIR_DISTANCE*vx;
			int j = _listeBases.get(i).getElementStructure();
			if (j>i)
			{
				drawLoop(i, j, x+(BASE_PAIR_DISTANCE*vx/2.0), y+(BASE_PAIR_DISTANCE*vy/2.0), dirAngle, coords, centers);
				centers[i].x = coords[i].x+BASE_PAIR_DISTANCE*vy;
				centers[i].y = y-BASE_PAIR_DISTANCE*vx;
				i = j;
				x += BASE_PAIR_DISTANCE*vx;
				y += BASE_PAIR_DISTANCE*vy;
				centers[i].x = coords[i].x+BASE_PAIR_DISTANCE*vy;
				centers[i].y = y-BASE_PAIR_DISTANCE*vx;
			}
			x += MULTILOOP_DISTANCE*vx;
			y += MULTILOOP_DISTANCE*vy;
			i += nextBase;
		}
	}
	
	
	/**
	 * Compute the symmetric of all the points in the points array
	 * relative to the line that goes through p and has director vector v.
	 * The array is modified in-place.
	 */
	private static void symmetric(
			Point2D.Double p,
			Point2D.Double v,
			Point2D.Double[] points) {
		// ||v||^2
		double lv = v.x*v.x + v.y*v.y;
		for (int i=0; i<points.length; i++) {
			// A is the coordinates of points[i] after moving the origin at p
			Point2D.Double A = new Point2D.Double(points[i].x - p.x, points[i].y - p.y);
			// Symmetric of A
			Point2D.Double B = new Point2D.Double(
					-(A.x*v.y*v.y - 2*A.y*v.x*v.y - A.x*v.x*v.x) / lv,
					 (A.y*v.y*v.y + 2*A.x*v.x*v.y - A.y*v.x*v.x) / lv);
			// Change the origin back
			points[i].x = B.x + p.x;
			points[i].y = B.y + p.y;
		}
	}
	
	/**
	 * Draw this RNA like the given template.
	 */
	public void drawRNATemplate(RNATemplate template) throws RNATemplateDrawingAlgorithmException {
		_drawn = true;
		_drawMode = DRAW_MODE_TEMPLATE;
		
		// debug
		try {
			RNA perfectMatchingRNA = template.toRNA();
			System.out.println("An RNA that would perfectly match this template would be:");
			System.out.println(perfectMatchingRNA.getStructDBN());
		} catch (ExceptionInvalidRNATemplate e) {
			e.printStackTrace();
		}
		
		RNATemplateMapping mapping = RNATemplateAlign.multiPassMap(this, template);
		
		// debug
//			RNATemplateAlign.printMapping(mapping, template, getSeq());
//			try {
//				TreeGraphviz.treeToGraphvizPostscript(alignment, "alignment_graphviz.ps");
//			} catch (IOException e) {
//				e.printStackTrace();
//			}
		
		// Allocate the coords and centers arrays
		// We create Point2D.Double objects in it but the algorithms
		// we use may choose to create new Point2D.Double objects or to
		// modify those created here.
		Point2D.Double[] coords = new Point2D.Double[_listeBases.size()];
		Point2D.Double[] centers = new Point2D.Double[_listeBases.size()];
		for (int i = 0; i < _listeBases.size(); i++) {
			coords[i] = new Point2D.Double(0, 0);
			centers[i] = new Point2D.Double(0, 0);
		}
		
		Set<RNATemplateHelix> alreadyDrawnHelixes = new HashSet<RNATemplateHelix>();
		RNATemplateHelix lastMappedHelix = null;
		EdgeEndPoint howWeGotOutOfLastHelix = null;
		int howWeGotOutOfLastHelixBaseIndex = -1;
		Iterator<RNATemplateElement> iter = template.rnaIterator();
		while (iter.hasNext()) {
			RNATemplateElement element = iter.next();
			if (element instanceof RNATemplateHelix
					&& mapping.getAncestor(element) != null) {
				// We have a mapping between an helix in the RNA sequence
				// and an helix in the template.
				
				RNATemplateHelix helix = (RNATemplateHelix) element;
				boolean firstTimeWeMeetThisHelix;
				int[] basesInHelixArray = RNATemplateAlign.intArrayFromList(mapping.getAncestor(helix));
				Arrays.sort(basesInHelixArray);
				
				// Draw this helix if it has not already been done
				if (!alreadyDrawnHelixes.contains(helix)) {
					firstTimeWeMeetThisHelix = true;
					drawHelixLikeTemplateHelix(basesInHelixArray, helix, coords, centers);
					alreadyDrawnHelixes.add(helix);
				} else {
					firstTimeWeMeetThisHelix = false;
				}
				
				EdgeEndPoint howWeGetInCurrentHelix;
				if (firstTimeWeMeetThisHelix) {
					if (helix.getIn1Is() == In1Is.IN1_IS_5PRIME) {
						howWeGetInCurrentHelix = helix.getIn1();
					} else {
						howWeGetInCurrentHelix = helix.getIn2();
					}
				} else {
					if (helix.getIn1Is() == In1Is.IN1_IS_5PRIME) {
						howWeGetInCurrentHelix = helix.getIn2();
					} else {
						howWeGetInCurrentHelix = helix.getIn1();
					}
				}
				
				if (lastMappedHelix != null) {
					// Now draw the RNA sequence (possibly containing helixes)
					// between the last template drawn helix and this one.
					Point2D.Double helixVector = new Point2D.Double();
					Point2D.Double P0 = new Point2D.Double();
					Point2D.Double P3 = new Point2D.Double();

					if (lastMappedHelix == helix) {
						// Last helix is the same as the current one so
						// nothing matched (or at best a single
						// non-paired sequence) so we will just
						// use the Radiate algorithm
						
						whereIsThisHelixEndpoint(howWeGotOutOfLastHelix, P0, helixVector, new Point2D.Double());
						whereIsThisHelixEndpoint(howWeGetInCurrentHelix, P3, new Point2D.Double(), new Point2D.Double());
						
						double angle = angleFromVector(helixVector);
						int b1 = basesInHelixArray[basesInHelixArray.length/2 - 1];
						int b2 = basesInHelixArray[basesInHelixArray.length/2];
						Point2D.Double loopCenter = new Point2D.Double((P0.x + P3.x)/2, (P0.y + P3.y)/2);
						drawLoop(b1,
								 b2,
								 loopCenter.x,
								 loopCenter.y,
								 angle,
								 coords,
								 centers);
						// If the helix is flipped, we need to compute the symmetric
						// of the whole loop.
						if (helix.isFlipped()) {
							Point2D.Double[] points1 = new Point2D.Double[b2-b1+1];
							Point2D.Double[] points2 = new Point2D.Double[b2-b1+1];
							for (int b=b1; b<=b2; b++) {
								points1[b-b1] = coords[b];
								points2[b-b1] = centers[b];
							}
							symmetric(loopCenter, helixVector, points1);
							symmetric(loopCenter, helixVector, points2);
						}
					} else {
						// No helixes matched between the last helix and
						// the current one, so we draw what is between
						// using the radiate algorithm but on the Bezier curve.
						
						Point2D.Double P1, P2;
						
						if (howWeGotOutOfLastHelix.getOtherElement() instanceof RNATemplateUnpairedSequence
								&& howWeGetInCurrentHelix.getOtherElement() instanceof RNATemplateUnpairedSequence) {
							// We will draw the bases on a Bezier curve
							P1 = new Point2D.Double();
							computeBezierPoints(howWeGotOutOfLastHelix, P0, P1, helixVector);
							
							P2 = new Point2D.Double();
							computeBezierPoints(howWeGetInCurrentHelix, P3, P2, new Point2D.Double());
						} else {
							// We will draw the bases on a straight line between P0 and P3
							P1 = null;
							P2 = null;
							whereIsThisHelixEndpoint(howWeGotOutOfLastHelix, P0, helixVector, new Point2D.Double());
							whereIsThisHelixEndpoint(howWeGetInCurrentHelix, P3, new Point2D.Double(), new Point2D.Double());
						}
						
						
						int b1 = howWeGotOutOfLastHelixBaseIndex;
						int b2;
						if (firstTimeWeMeetThisHelix) {
							b2 = basesInHelixArray[0];
						} else {
							b2 = basesInHelixArray[basesInHelixArray.length/2];
						}
						drawAlongBezierCurve(b1+1, b2-1, P0, P1, P2, P3, coords, centers);
					}
					
				} else {
					// here we draw what is before the first mapped helix
					//drawFromHelixEndPoint(howWeGetInCurrentHelix, basesInHelixArray[0], coords, centers);
				}
				
				lastMappedHelix = helix;
				howWeGotOutOfLastHelix = howWeGetInCurrentHelix.getNextEndPoint();
				if (firstTimeWeMeetThisHelix) {
					howWeGotOutOfLastHelixBaseIndex = basesInHelixArray[basesInHelixArray.length/2-1];
				} else {
					howWeGotOutOfLastHelixBaseIndex = basesInHelixArray[basesInHelixArray.length-1];
				}
			}
		} // end template iteration
		
		// now we should draw what remains after the last mapped helix
		// TODO
		
		// Now move all bases, according to arrays coords and centers
		// and taking in account the space between bases parameter.
		for (int i = 0; i < _listeBases.size(); i++) {
			_listeBases.get(i).setCoords(
					new Point2D.Double(coords[i].x * _spaceBetweenBases,
							coords[i].y * _spaceBetweenBases));
			_listeBases.get(i).setCenter(
					new Point2D.Double(centers[i].x * _spaceBetweenBases,
							centers[i].y * _spaceBetweenBases));
		}
	}
	
	



	private double objFun(int n1, int n2, double r) {
		return (((double) n1) * 2.0
				* Math.asin(((double) BASE_PAIR_DISTANCE) / (2.0 * r))
				+ ((double) n2) * 2.0
				* Math.asin(((double) MULTILOOP_DISTANCE) / (2.0 * r)) - (2.0 * Math.PI));
	}

	private double determineRadius(int n1, int n2, double startRadius) {
		double xmin = BASE_PAIR_DISTANCE / 2.0;
		double xmax = 3.0 * BASE_PAIR_DISTANCE + 1;
		double x = (xmin + xmax) / 2.0;
		double y = 10000.0;
		double ymin = -1000.0;
		double ymax = 1000.0;
		int numIt = 0;
		double precision = 0.00001;
		while ((Math.abs(y) > precision) && (numIt < 10000)) {
			x = (xmin + xmax) / 2.0;
			y = objFun(n1, n2, x);
			ymin = objFun(n1, n2, xmax);
			ymax = objFun(n1, n2, xmin);
			if (ymin > 0.0) {
				xmax = xmax + (xmax - xmin);
			} else if ((y <= 0.0) && (ymax > 0.0)) {
				xmax = x;
			} else if ((y >= 0.0) && (ymin < 0.0)) {
				xmin = x;
			} else if (ymax < 0.0) {
				xmin = Math.max(xmin - (x - xmin), Math.max(
						BASE_PAIR_DISTANCE / 2.0, MULTILOOP_DISTANCE / 2.0));
				xmax = x;
			}
			numIt++;
		}
		// if (Math.abs(y)>precision)
		// System.out.println("x: ["+xmin+","+x+","+xmax+"] inc:"+increment+" y: ["+ymin+","+y+","+ymax+"]");
		return x;
	}

	public void drawRNA() throws ExceptionNAViewAlgorithm {
		drawRNA(RNA.DEFAULT_DRAW_MODE);
	}
	public void drawRNA(int mode) throws ExceptionNAViewAlgorithm {
		_drawMode = mode;
		switch (get_drawMode()) {
		case RNA.DRAW_MODE_RADIATE:
			drawRNARadiate();
			break;
		case RNA.DRAW_MODE_LINEAR:
			drawRNALine();
			break;
		case RNA.DRAW_MODE_CIRCULAR:
			drawRNACircle();
			break;
		case RNA.DRAW_MODE_NAVIEW:
			drawRNANAView();
			break;
		case RNA.DRAW_MODE_VARNA_VIEW:
			drawRNAVARNAView();
			break;
		case RNA.DRAW_MODE_MOTIFVIEW:
			drawMOTIFView();
			break;
		default:
			break;
		}

	}
	

	public int getDrawMode() {
		return _drawMode;
	}

	
	
	private void drawLoop(int i, int j, double x, double y, double dirAngle,
			Point2D.Double[] coords, Point2D.Double[] centers) {
		if (i > j) {
			return;
		}

		// BasePaired
		if (_listeBases.get(i).getElementStructure() == j) {
			double normalAngle = Math.PI / 2.0;
			centers[i] = new Point2D.Double(x, y);
			centers[j] = new Point2D.Double(x, y);
			coords[i].x = (x + BASE_PAIR_DISTANCE
					* Math.cos(dirAngle - normalAngle) / 2.0);
			coords[i].y = (y + BASE_PAIR_DISTANCE
					* Math.sin(dirAngle - normalAngle) / 2.0);
			coords[j].x = (x + BASE_PAIR_DISTANCE
					* Math.cos(dirAngle + normalAngle) / 2.0);
			coords[j].y = (y + BASE_PAIR_DISTANCE
					* Math.sin(dirAngle + normalAngle) / 2.0);
			drawLoop(i + 1, j - 1, x + LOOP_DISTANCE * Math.cos(dirAngle), y
					+ LOOP_DISTANCE * Math.sin(dirAngle), dirAngle, coords,
					centers);
		} else {
			int k = i;
			Vector<Integer> basesMultiLoop = new Vector<Integer>();
			Vector<Integer> helices = new Vector<Integer>();
			int l;
			while (k <= j) {
				l = _listeBases.get(k).getElementStructure();
				if (l > k) {
					basesMultiLoop.add(new Integer(k));
					basesMultiLoop.add(new Integer(l));
					helices.add(new Integer(k));
					k = l + 1;
				} else {
					basesMultiLoop.add(new Integer(k));
					k++;
				}
			}
			int mlSize = basesMultiLoop.size() + 2;
			int numHelices = helices.size() + 1;
			double totalLength = MULTILOOP_DISTANCE * (mlSize - numHelices)
					+ BASE_PAIR_DISTANCE * numHelices;
			double multiLoopRadius;
			double angleIncrementML;
			double angleIncrementBP;
			if (mlSize > 3) {
				multiLoopRadius = determineRadius(numHelices, mlSize
						- numHelices, (totalLength) / (2.0 * Math.PI));
				angleIncrementML = -2.0
						* Math.asin(((float) MULTILOOP_DISTANCE)
								/ (2.0 * multiLoopRadius));
				angleIncrementBP = -2.0
						* Math.asin(((float) BASE_PAIR_DISTANCE)
								/ (2.0 * multiLoopRadius));
			} else {
				multiLoopRadius = 35.0;
				angleIncrementBP = -2.0
						* Math.asin(((float) BASE_PAIR_DISTANCE)
								/ (2.0 * multiLoopRadius));
				angleIncrementML = (-2.0 * Math.PI - angleIncrementBP) / 2.0;
			}
			// System.out.println("MLr:"+multiLoopRadius+" iBP:"+angleIncrementBP+" iML:"+angleIncrementML);

			double centerDist = Math.sqrt(Math.max(Math.pow(multiLoopRadius, 2)
					- Math.pow(BASE_PAIR_DISTANCE / 2.0, 2), 0.0))
					- LOOP_DISTANCE;
			Point2D.Double mlCenter = new Point2D.Double(
					(x + (centerDist * Math.cos(dirAngle))),
					(y + (centerDist * Math.sin(dirAngle))));

			// Base directing angle for (multi|hairpin) loop, from the center's
			// perspective
			double baseAngle = dirAngle
			// U-turn
					+ Math.PI
					// Account for already drawn supporting base-pair
					+ 0.5 * angleIncrementBP
					// Base cannot be paired twice, so next base is at
					// "unpaired base distance"
					+ 1.0 * angleIncrementML;
			double[] angles = new double[_listeBases.size()];
			int n1 = 1;
			int n2 = 1;
			for (k = basesMultiLoop.size() - 1; k >= 0; k--) {
				l = basesMultiLoop.get(k).intValue();
				centers[l] = mlCenter;
				angles[l] = baseAngle;
				coords[l].x = mlCenter.x + multiLoopRadius
						* Math.cos(baseAngle);
				coords[l].y = mlCenter.y + multiLoopRadius
						* Math.sin(baseAngle);
				if ((_listeBases.get(l).getElementStructure() < l)
						&& (_listeBases.get(l).getElementStructure() != -1)) {
					baseAngle += angleIncrementBP;
					n1++;
				} else {
					baseAngle += angleIncrementML;
					n2++;
				}
			}
			// System.out.println("n1:"+n1+" n2:"+n2);
			double newAngle;
			int m, n;
			for (k = 0; k < helices.size(); k++) {
				m = helices.get(k).intValue();
				n = _listeBases.get(m).getElementStructure();
				newAngle = (angles[m] + angles[n]) / 2.0;
				drawLoop(m + 1, n - 1, (LOOP_DISTANCE * Math.cos(newAngle))
						+ (coords[m].x + coords[n].x) / 2.0,
						(LOOP_DISTANCE * Math.sin(newAngle))
								+ (coords[m].y + coords[n].y) / 2.0, newAngle,
						coords, centers);
			}
		}
	}


	public void drawRNARadiate() {
		drawRNARadiate(-1.0,_flatExteriorLoop);
	}
	
	public void drawRNARadiate(double dirAngle, boolean flatExteriorLoop) {
		_drawn = true;
		_drawMode = DRAW_MODE_RADIATE;
		Point2D.Double[] coords = new Point2D.Double[_listeBases.size()];
		Point2D.Double[] centers = new Point2D.Double[_listeBases.size()];
		for (int i = 0; i < _listeBases.size(); i++) {
			coords[i] = new Point2D.Double(0, 0);
			centers[i] = new Point2D.Double(0, 0);
		}
		if (flatExteriorLoop)
		{
		  dirAngle += 1.0 - Math.PI/2.0;
		  int i=0;
		  double x = 0.0;
		  double y = 0.0;
		  double vx = -Math.sin(dirAngle);
		  double vy = Math.cos(dirAngle);
		  while(i<_listeBases.size())
		  {
			  coords[i].x = x; 
			  coords[i].y = y;
			  centers[i].x = x+BASE_PAIR_DISTANCE*vy;
			  centers[i].y = y-BASE_PAIR_DISTANCE*vx;
			  int j = _listeBases.get(i).getElementStructure();
			  if (j>i)
			  {
				  drawLoop(i, j, x+(BASE_PAIR_DISTANCE*vx/2.0), y+(BASE_PAIR_DISTANCE*vy/2.0), dirAngle, coords, centers);
				  centers[i].x = coords[i].x+BASE_PAIR_DISTANCE*vy;
				  centers[i].y = y-BASE_PAIR_DISTANCE*vx;
				  i = j;
				  x += BASE_PAIR_DISTANCE*vx;
				  y += BASE_PAIR_DISTANCE*vy;
				  centers[i].x = coords[i].x+BASE_PAIR_DISTANCE*vy;
				  centers[i].y = y-BASE_PAIR_DISTANCE*vx;
			  }
			  x += MULTILOOP_DISTANCE*vx;
			  y += MULTILOOP_DISTANCE*vy;
			  i += 1;
		  }
		}
		else
		{
		  drawLoop(0, _listeBases.size() - 1, 0, 0, dirAngle, coords, centers);
		}
		for (int i = 0; i < _listeBases.size(); i++) {
			_listeBases.get(i).setCoords(
					new Point2D.Double(coords[i].x * _spaceBetweenBases,
							coords[i].y * _spaceBetweenBases));
			_listeBases.get(i).setCenter(
					new Point2D.Double(centers[i].x * _spaceBetweenBases,
							centers[i].y * _spaceBetweenBases));
		}

		// TODO
		// change les centres des bases de la premiere helice vers la boucle la
		// plus proche
	}

	public void drawRNANAView() throws ExceptionNAViewAlgorithm {
		_drawMode = DRAW_MODE_NAVIEW;
		_drawn = true;

		ArrayList<Double> X = new ArrayList<Double>(_listeBases.size());
		ArrayList<Double> Y = new ArrayList<Double>(_listeBases.size());
		ArrayList<Short> pair_table = new ArrayList<Short>(_listeBases.size());

		for (int i = 0; i < _listeBases.size(); i++) {
			pair_table.add(Short.valueOf(String.valueOf(_listeBases.get(i)
					.getElementStructure())));
		}
		NAView naView = new NAView();
		naView.naview_xy_coordinates(pair_table, X, Y);

		// Updating individual base positions
		for (int i = 0; i < _listeBases.size(); i++) {
			_listeBases.get(i).setCoords(
					new Point2D.Double(X.get(i) * 2.5 * _spaceBetweenBases, Y
							.get(i)
							* 2.5 * _spaceBetweenBases));
		}

		// Updating centers
		for (int i = 0; i < _listeBases.size(); i++) {
			int indicePartner = _listeBases.get(i).getElementStructure();
			if (indicePartner != -1) {
				Point2D.Double base = _listeBases.get(i).getCoords();
				Point2D.Double partner = _listeBases.get(indicePartner)
						.getCoords();
				_listeBases.get(i).setCenter(
						new Point2D.Double((base.x + partner.x) / 2.0,
								(base.y + partner.y) / 2.0));
			} 
			else {
				Vector<Integer> loop = getLoopBases(i);
				double tmpx = 0.0;
				double tmpy = 0.0;
				for (int j = 0; j < loop.size(); j++) {
					int partner = loop.elementAt(j);
					Point2D.Double loopmember = _listeBases.get(partner)
							.getCoords();
					tmpx += loopmember.x;
					tmpy += loopmember.y;
				}
				_listeBases.get(i).setCenter(
						new Point2D.Double(tmpx / loop.size(), tmpy
								/ loop.size()));
			}
		}
	}

	public void drawMOTIFView() {
		_drawn = true;
		_drawMode = DRAW_MODE_MOTIFVIEW;	
		int spaceBetweenStrand =0;
		Motif motif = new Motif(this,get_listeBases());
		motif.listStrand();
		for (int i = 0; i < motif.getListStrand().sizeStruct(); i++ ){
			for (int j = 0; j < motif.getListStrand().getStrand(i).sizeStrand(); j++ ){
				int indice = motif.getListStrand().getStrand(i).getMB(j).getIndex();
				get_listeBases().get(indice).setCoords(
						new Point2D.Double(0,0));
				get_listeBases().get(indice).setCenter(
						new Point2D.Double(0, 0));

			}
		}
		//Recherche du brin central
		int centralStrand = motif.getCentralStrand();
		
		//Cas o� l'on a un motif en �toile
		if(centralStrand!=-1){
			//On positionne le brin central
			motif.positionneSpecificStrand(centralStrand, spaceBetweenStrand);

			//On place les autres brins par rapport a ce brin central
			motif.orderStrands(centralStrand);
		}
		
		else {	
			centralStrand = 0;
			motif.positionneStrand();
			motif.ajusteStrand();
		}
		motif.reajustement();
		motif.deviationBasePair();	
		motif.setCenterMotif();
	}
	
		
	public ArrayList<ModeleBase> getAllPartners(int indice) {		
		ArrayList<ModeleBase> result = new ArrayList<ModeleBase>();
		ModeleBase me = this.getBaseAt(indice);
		int i = me.getElementStructure();
		if (i!= -1)
		{ result.add(getBaseAt(i)); }
		ArrayList<ModeleStyleBP> msbps = getAuxBPs(indice);
		for (ModeleStyleBP m : msbps)
		{ result.add(m.getPartner(me)); }		
		return result;
	}
	
	
	public int get_drawMode() {
		return _drawMode;
	}

	public void setDrawMode(int drawMode) {
		_drawMode = drawMode;
	}

	public void setRNA(char[] seq, int[] str)
		throws ExceptionFileFormatOrSyntax 
	{
		setRNA(seq, str, 1);
	}
	
	public void setRNA(char[] seq, int[] str, int baseIndex)
			throws ExceptionFileFormatOrSyntax {
		clearAnnotations();
		_listeBases = new ArrayList<ModeleBase>();
		if (seq.length != str.length) {
			warningEmition("Sequence length " + seq.length
					+ " differs from that of secondary structure " + str.length
					+ ". \nAdapting sequence length ...");
			if (seq.length < str.length) {
				String seqtmp = String.copyValueOf(seq);
				while (seqtmp.length() < str.length) {
					seqtmp += " ";
				}
				seq = seqtmp.toCharArray();
			} else {
				seq = String.copyValueOf(seq).substring(0, str.length)
						.toCharArray();
			}
		}
		for (int i = 0; i < str.length; i++) {
			_listeBases.add(new ModeleBaseNucleotide(seq[i], i, baseIndex+i));
		}
		applyStruct(str);
	}

	/**
	 * Sets the RNA to be drawn. Uses when comparison mode is on. Will draw the
	 * super-structure passed in parameters and apply specials styles to the
	 * bases owning by each RNA alignment and both.
	 * 
	 * @param seq
	 *            - The sequence of the super-structure This sequence shall be
	 *            designed like this:
	 *            <code>firstRNA1stBaseSecondRNA1stBaseFirstRNA2ndBaseSecondRNA2ndBase [...]</code>
	 * <br>
	 *            <b>Example:</b> <code>AAC-GUAGA--UGG</code>
	 * @param struct
	 *            - The super-structure
	 * @param basesOwn
	 *            - The RNA owning bases array (each index will be:0 when common
	 *            base, 1 when first RNA alignment base, 2 when second RNA
	 *            alignment base)
	 * @throws ExceptionUnmatchedClosingParentheses
	 * @throws ExceptionFileFormatOrSyntax 
	 */
	public void setRNA(String seq, String struct, ArrayList<Integer> basesOwn)
			throws ExceptionUnmatchedClosingParentheses, ExceptionFileFormatOrSyntax {
		clearAnnotations();
		_listeBases = new ArrayList<ModeleBase>();
		// On "parse" la structure (repérage des points, tiret et couples
		// parentheses ouvrante/fermante)
		int[] array_struct = parseStruct(struct);
		// On enregistre la taille de la structure pour éviter d'appeler la
		// fonction à chaque passage dans la boucle for ci-après
		int size = struct.length();
		// Compteur pour les caractères (on va aller de deux en deux)
		int j = 0;
		boolean simple = false;
		for (int i = 0; i < size; i++) {
			// On créer une nouvelle base de comparaison avec les deux bases
			// courantes
			// if (seq.charAt(j) != seq.charAt(j + 1)) {
			_listeBases.add(new ModeleBasesComparison(seq.charAt(j), seq
					.charAt(j + 1), i));
			// }
			// Si les deux bases sont les mêmes, alors on dessinne une base
			// simple
			// else {
			// _listeBases.add(new ModeleBaseNucleotide(seq.charAt(j)));
			// simple = true;
			// }
			// On applique la structure à cette base
			_listeBases.get(i).setElementStructure(array_struct[i]);
			// Non-testé: On renseigne à quelle ARN cette base appartient à
			// l'origine (en termes de structure).
			if (!simple)
				((ModeleBasesComparison) _listeBases.get(i))
						.set_appartenance(basesOwn.get(i));
			// S'il s'agit d'une base simple
			// else {
			// System.out.println("This is "+_listeBases.get(i)+" for "+i);
			// if (basesOwn.get(i) == 0) {
			// ((ModeleBaseNucleotide)
			// _listeBases.get(i)).get_styleBase().set_base_inner_color(ModeleBasesComparison.BOTH_RNA_COLOR);
			// }
			// else if (basesOwn.get(i) == 1){
			// ((ModeleBaseNucleotide)
			// _listeBases.get(i)).get_styleBase().set_base_inner_color(ModeleBasesComparison.FIRST_RNA_COLOR);
			// }
			// else if ((basesOwn.get(i) == 2)) {
			// ((ModeleBaseNucleotide)
			// _listeBases.get(i)).get_styleBase().set_base_inner_color(ModeleBasesComparison.SECOND_RNA_COLOR);
			// }
			// }
			// On passe au morceau de la séquence suivante.
			j += 2;
		}
	}

	public void setRNA(String seq, int[] str)
			throws ExceptionFileFormatOrSyntax {
		setRNA(seq.toCharArray(), str);
	}

	public void setRNA(String seq, String dbnStr)
			throws ExceptionUnmatchedClosingParentheses,
			ExceptionFileFormatOrSyntax {
		clearAnnotations();
		String parDBN = dbnStr.replace('(', '(').replace(')', ')').replace('[', ':').replace(']', ':').replace('{', ':').replace('}', ':');
		String braDBN = dbnStr.replace('(', ':').replace(')', ':').replace('[', '(').replace(']', ')').replace('{', ':').replace('}', ':');
		String accDBN = dbnStr.replace('(', ':').replace(')', ':').replace('[', ':').replace(']', ':').replace('{', '(').replace('}', ')');
	    int[] parStr = parseStruct(parDBN);
		int[] braStr = parseStruct(braDBN);
		int[] accStr = parseStruct(accDBN);
		int[] finStr = new int[parStr.length];
		for(int i=0;i<parStr.length;i++)
		{ finStr[i] = -1; } 
		
		for(int i=0;i<parStr.length;i++)
		{
			if (parStr[i]>i)
			{
				finStr[i] = parStr[i];
				finStr[finStr[i]] = i;
			}
			else if (braStr[i]>i)
			{
				if ((parStr[i]==-1)&&(parStr[braStr[i]]==-1))
				{
					finStr[i] = braStr[i];
					finStr[finStr[i]] = i;				
				}
			}
			else if (accStr[i]>i)
			{
				if ((parStr[i]==-1)&&(parStr[accStr[i]]==-1)&&(braStr[i]==-1)&&(braStr[accStr[i]]==-1))
				{
					finStr[i] = accStr[i];
					finStr[finStr[i]] = i;				
				}
			}
		}
		setRNA(seq.toCharArray(), finStr);
	}

	public int[] parseStruct(String str)
			throws ExceptionUnmatchedClosingParentheses, ExceptionFileFormatOrSyntax {
		int[] result = new int[str.length()];
		int unexpectedChar = -1;
		Stack<Integer> p = new Stack<Integer>();
		for (int i = 0; i < str.length(); i++) {
			char c = str.charAt(i);
			if (c == '(') {
				p.push(new Integer(i));
			} else if (c == '.' || c == '-' || c == ':') {
				result[i] = -1;
			} else if (c == ')') {
				if (p.size() == 0) {
					throw new ExceptionUnmatchedClosingParentheses(i + 1);
				}
				int j = p.pop().intValue();
				result[i] = j;
				result[j] = i;
			} else {
				result[i] = -1;
				if (unexpectedChar == -1)
					unexpectedChar = i;
			}
		}

		if (unexpectedChar != -1) {
			//warningEmition("Unexpected Character at index:" + unexpectedChar);
			throw new ExceptionFileFormatOrSyntax("Unexpected Character at index:" + unexpectedChar); 
		}

		if (p.size() != 0) {
			throw new ExceptionUnmatchedClosingParentheses(
					p.pop().intValue() + 1);
		}

		return result;
	}

	public Point getHelixInterval(int index) {
		if ((index < 0) || (index >= _listeBases.size())) {
			return new Point(index, index);
		}
		int j = _listeBases.get(index).getElementStructure();
		if (j != -1) {
			int minH = index;
			int maxH = index;
			if (j > index) {
				maxH = j;
			} else {
				minH = j;
			}
			boolean over = false;
			while (!over) {
				if ((minH < 0) || (maxH >= _listeBases.size())) {
					over = true;
				} else {
					if (_listeBases.get(minH).getElementStructure() == maxH) {
						minH--;
						maxH++;
					} else {
						over = true;
					}
				}
			}
			minH++;
			maxH--;
			return new Point(minH, maxH);
		}
		return new Point(0, 0);
	}

	public ArrayList<Integer> getHelix(int index) {
		ArrayList<Integer> result = new ArrayList<Integer>();
		if ((index < 0) || (index >= _listeBases.size())) {
			return result;
		}
		Point p = getHelixInterval(index);
		for (int i =p.x;i<=p.y;i++)
		{
			result.add(i);
			result.add(this._listeBases.get(i).getElementStructure());
		}
		return result;
	}

	public Point getMultiLoop(int index) {
		if ((index < 0) || (index >= _listeBases.size())) {
			return new Point(index, index);
		}
		Point h = getHelixInterval(index);
		int minH = h.x - 1;
		int maxH = h.y + 1;
		boolean over = false;
		while (!over) {
			if (minH < 0) {
				over = true;
				minH = 0;
			} else {
				if (_listeBases.get(minH).getElementStructure() == -1) {
					minH--;
				} else if (_listeBases.get(minH).getElementStructure() < minH) {
					minH = _listeBases.get(minH).getElementStructure() - 1;
				} else {
					over = true;
				}
			}
		}
		over = false;
		while (!over) {
			if (maxH > _listeBases.size() - 1) {
				over = true;
				maxH = _listeBases.size() - 1;
			} else {
				if (_listeBases.get(maxH).getElementStructure() == -1) {
					maxH++;
				} else if (_listeBases.get(maxH).getElementStructure() > maxH) {
					maxH = _listeBases.get(maxH).getElementStructure() + 1;
				} else {
					over = true;
				}
			}
		}
		return new Point(minH, maxH);
	}

	public Vector<Integer> getLoopBases(int startIndex) {
		Vector<Integer> result = new Vector<Integer>();

		if ((startIndex < 0) || (startIndex >= _listeBases.size())) {
			return result;
		}
		int index = startIndex;
		result.add(startIndex);
		if (_listeBases.get(index).getElementStructure() <= index) {
			index = (index + 1) % _listeBases.size();
		} else {
			index = _listeBases.get(index).getElementStructure();
			result.add(index);
			index = (index + 1) % _listeBases.size();
		}

		while (index != startIndex) {
			result.add(index);
			if (_listeBases.get(index).getElementStructure() == -1) {
				index = (index + 1) % _listeBases.size();
			} else {
				index = _listeBases.get(index).getElementStructure();
				result.add(index);
				index = (index + 1) % _listeBases.size();
			}
		}
		return result;
	}

	/**
	 * Returns the RNA secondary structure displayed by this panel as a
	 * well-parenthesized word, accordingly to the DBN format
	 * 
	 * @return This panel's secondary structure
	 */
	public String getStructDBN() {
		String result = "";
		for (int i = 0; i < _listeBases.size(); i++) {
			int j = _listeBases.get(i).getElementStructure();
			if (j == -1) {
				result += ".";
			} else if (i > j) {
				result += ")";
			} else {
				result += "(";
			}
		}
		return result;
	}

	public String getStructDBN(int[] str) {
		String result = "";
		for (int i = 0; i < str.length; i++) {
			if (str[i] == -1) {
				result += ".";
			} else if (str[i] > i) {
				result += "(";
			} else {
				result += ")";
			}
		}
		return result;
	}

	/**
	 * Returns the raw nucleotides sequence for the displayed RNA
	 * 
	 * @return The RNA sequence
	 */
	public String getSeq() {
		String result = "";
		if (_comparisonMode) {
			for (int i = 0; i < _listeBases.size(); i++) {
				result += ((ModeleBasesComparison) _listeBases.get(i))
						.getBases();
			}
		} else {
			for (int i = 0; i < _listeBases.size(); i++) {
				result += ((ModeleBaseNucleotide) _listeBases.get(i)).get_c();
			}
		}
		return result;
	}

	public String getStructBPSEQ() {
		String result = "";
		if (_comparisonMode) {
			for (int i = 0; i < _listeBases.size(); i++) {
				result += (i + 1)
						+ " "
						+ ((ModeleBasesComparison) _listeBases.get(i))
								.getBases() + " "
						+ (_listeBases.get(i).getElementStructure() + 1)
						+ "\n";
			}
		} else {
			int[] str = getNonOverlappingStruct();
			for (int i = 0; i < _listeBases.size(); i++) {
				result += (i + 1) + " "
						+ ((ModeleBaseNucleotide) _listeBases.get(i)).get_c()
						+ " " + (str[i] + 1)
						+ "\n";
			}
		}
		return result;
	}
	
	public int[] getNonCrossingStruct()
	{
	  int[] result = new int[_listeBases.size()];
      // Adding "planar" base-pairs
	  for (int i = 0; i < _listeBases.size(); i++) {
			result[i] = _listeBases.get(i).getElementStructure();
		}
	  return result;
	}

	public int[] getNonOverlappingStruct()
	{
	  int[] result = getNonCrossingStruct();
	  // Adding additional base pairs when possible (No more than one base-pair per base)
	  for (int i = 0; i < _structureAux.size(); i++) {
			ModeleStyleBP msbp = _structureAux.get(i);
			ModeleBase mb5 = msbp.getPartner5();
			ModeleBase mb3 = msbp.getPartner3();
			int j5 = mb5.getIndex();
			int j3 = mb3.getIndex();
			if ((result[j3]==-1) && (result[j5]==-1))
			{
				result[j3] = j5;
				result[j5] = j3;
			}
		}
	  return result;
	}
	
	
	public String getStructCT() {
		String result = "";
		if (_comparisonMode) {
			for (int i = 0; i < _listeBases.size(); i++) {
				result += (i + 1)
						+ " "
						+ ((ModeleBasesComparison) _listeBases.get(i))
								.getBases() + " " + i + " " + (i + 2) + " "
						+ (_listeBases.get(i).getElementStructure() + 1) + " "
						+ (i + 1) + "\n";
			}
		} else {
			int[] str = getNonOverlappingStruct();
			for (int i = 0; i < _listeBases.size(); i++) {
				result += (i + 1) + " "
						+ ((ModeleBaseNucleotide) _listeBases.get(i)).get_c()
						+ " " + i + " " + (i + 2) + " "
						+ (str[i] + 1) + " "
						+ (i + 1) + "\n";
			}
		}
		return result;
	}

	public void saveAsBPSEQ(String path, String title)
			throws ExceptionExportFailed, ExceptionPermissionDenied {
		try {
			FileWriter f = new FileWriter(path);
			f.write("# " + title + "\n");
			f.write(this.getStructBPSEQ() + "\n");
			f.close();
		} catch (IOException e) {
			throw new ExceptionExportFailed(e.getMessage(), path);
		}
	}

	public void saveAsCT(String path, String title)
			throws ExceptionExportFailed, ExceptionPermissionDenied {
		try {
			FileWriter f = new FileWriter(path);
			f.write("" + _listeBases.size() + " " + title + "\n");
			f.write(this.getStructCT() + "\n");
			f.close();
		} catch (IOException e) {
			throw new ExceptionExportFailed(e.getMessage(), path);
		}
	}

	public void saveAsDBN(String path, String title)
			throws ExceptionExportFailed, ExceptionPermissionDenied {
		try {
			FileWriter f = new FileWriter(path);
			f.write("> " + title + "\n");
			f.write(getListeBasesToString() + "\n");
			f.write(getStructDBN() + "\n");
			f.close();
		} catch (IOException e) {
			throw new ExceptionExportFailed(e.getMessage(), path);
		}
	}

	public String getListeBasesToString() {
		String s = new String();
		if (_comparisonMode) {
			for (int i = 0; i < _listeBases.size(); i++) {
				s += ((ModeleBasesComparison) _listeBases.get(i)).getBases();
			}
		} else {
			for (int i = 0; i < _listeBases.size(); i++) {
				s += ((ModeleBaseNucleotide) _listeBases.get(i)).get_c();
			}
		}
		return s;
	}

	private boolean loadSecStrBPSEQ(Reader r) throws ExceptionPermissionDenied,
			ExceptionLoadingFailed, ExceptionFileFormatOrSyntax {
		boolean loadOk = false;
		try {
			BufferedReader fr = new BufferedReader(r);
			String line = fr.readLine();
			String seqTmp = "";
			Vector<Integer> strTmp = new Vector<Integer>();

			int bpFrom;
			char base;
			int bpTo;
			int minIndex = -1;
			boolean noWarningYet = true;
			String title = "";
			String filenameStr = "Filename:";
			String organismStr = "Organism:";
			String ANStr = "Accession Number:";
			while (line != null) {
				line = line.trim();
				String[] tokens = line.split("\\s+");
				if (tokens.length == 3 && !tokens[0].contains("#")&& !line.startsWith("Organism:")
						 && !line.startsWith("Filename:")&& !line.startsWith("Accession Number:")) 
				{
						base = tokens[1].charAt(tokens[1].length() - 1);

						bpFrom = (Integer.parseInt(tokens[0]));
						bpTo = (Integer.parseInt(tokens[2]));
						
						if (minIndex<0)
    						minIndex = bpFrom;
						bpFrom -= minIndex; 
						if (bpTo!=0)
						  bpTo -= minIndex;
						else
						  bpTo = -1;

						if (bpFrom != seqTmp.length()) {
							if (noWarningYet) {
								noWarningYet = false;
								warningEmition("Discontinuity detected between nucleotides "
										+ (seqTmp.length())
										+ " and "
										+ (bpFrom + 1)
										+ "!\nFilling in missing portions with unpaired unknown 'X' nucleotides ...");
							}
							while (bpFrom != seqTmp.length()) {
								seqTmp += 'X';
								strTmp.add(-1);
							}
						}
						seqTmp += base;
						strTmp.add(bpTo);
				}
				else if (tokens[0].startsWith("#"))
				{
					int occur = line.indexOf("#");
					String tmp = line.substring(occur+1);
					title += tmp.trim()+" ";
				}
				else if (tokens[0].startsWith(filenameStr))
				{
					int occur = line.indexOf(filenameStr);
					String tmp = line.substring(occur+filenameStr.length());
					title += tmp.trim();
				}
				else if (tokens[0].startsWith(organismStr))
				{
					int occur = line.indexOf(organismStr);
					String tmp = line.substring(occur+organismStr.length());
					if (title.length()!=0)
					{
						title = "/"+title;
					}
					title = tmp.trim() + title;
				}
				else if (line.contains(ANStr))
				{
					int occur = line.indexOf(ANStr);
					String tmp = line.substring(occur+ANStr.length());
					if (title.length()!=0)
					{
						title += " ";
					}
					title +="("+tmp.trim()+")";
				}
				line = fr.readLine();
			}
			if (strTmp.size() != 0) {
				char[] seq = seqTmp.toCharArray();
				int[] str = new int[strTmp.size()];
				for (int i = 0; i < strTmp.size(); i++) {
					str[i] = strTmp.elementAt(i).intValue();
				}
				
				setRNA(seq, str, minIndex);
				setName(title);
				loadOk = true;
			}
		}
		catch (NumberFormatException e) {
			return false;
		}
		catch (Exception e) {
			throw new ExceptionLoadingFailed(e.getMessage(), "");
		} 
		return loadOk;
	}

	private boolean loadSecStrCT(Reader r) throws ExceptionPermissionDenied,
			ExceptionLoadingFailed, ExceptionFileFormatOrSyntax {
		boolean loadOk = false;
		try {
			BufferedReader fr = new BufferedReader(r);
			String line = fr.readLine();
			String seqTmp = "";
			Vector<Integer> strTmp = new Vector<Integer>();
			int bpFrom;
			char base;
			int bpTo;
			boolean noWarningYet = true;
			int minIndex = -1;
			String title = "";
			while (line != null) {
				line = line.trim();
				String[] tokens = line.split("\\s+");
				if (tokens.length >= 6) {
					try{
					bpFrom = (Integer.parseInt(tokens[0]));
					bpTo = (Integer.parseInt(tokens[4]));
					if (minIndex==-1)
						minIndex = bpFrom;
					bpFrom -= minIndex;
					if (bpTo!=0)
						bpTo -= minIndex;
					else
						bpTo = -1;
					base = tokens[1].charAt(tokens[1].length() - 1);
					Integer.parseInt(tokens[2]);
					Integer.parseInt(tokens[3]);
					Integer.parseInt(tokens[5]);
					if (bpFrom != seqTmp.length()) {
						if (noWarningYet) {
							noWarningYet = false;
							warningEmition("Discontinuity detected between nucleotides "
									+ (seqTmp.length())
									+ " and "
									+ (bpFrom + 1)
									+ "!\nFilling in missing portions with unpaired unknown 'X' nucleotides ...");
						}
						while (bpFrom != seqTmp.length()) {
							seqTmp += 'X';
							strTmp.add(-1);
						}
					}
					seqTmp += base;
					strTmp.add(bpTo);
					}
					catch (NumberFormatException e) {
						}
				}
				if ((line.contains("ENERGY = "))||line.contains("dG = ")) 
				{
					String[] ntokens = line.split("\\s+");
					if (ntokens.length>=4)
					{
						String energy = ntokens[3];
						for(int i=4;i<ntokens.length;i++)
						{
							title += ntokens[i]+" ";
						}
						title += "(E="+energy+" KCal/Mol)";
					}
				}
				line = fr.readLine();
			}
			if (strTmp.size() != 0) {
				char[] seq = seqTmp.toCharArray();
				int[] str = new int[strTmp.size()];
				for (int i = 0; i < strTmp.size(); i++) {
					str[i] = strTmp.elementAt(i).intValue();
				}
				setRNA(seq, str, minIndex);
				setName(title);
				loadOk = true;
			}
		} catch (IOException e) {
			throw new ExceptionLoadingFailed(e.getMessage(), "");
		} catch (NumberFormatException e) {
			throw new ExceptionFileFormatOrSyntax(e.getMessage(), "");
		}
		return loadOk;
	}

	private boolean loadSecStrRNAML(Reader r) throws ExceptionPermissionDenied,
			ExceptionLoadingFailed, ExceptionFileFormatOrSyntax {
		boolean loadOk = false;
		try {
			System.setProperty("javax.xml.parsers.SAXParserFactory", "com.sun.org.apache.xerces.internal.jaxp.SAXParserFactoryImpl");
			SAXParserFactory saxFact = javax.xml.parsers.SAXParserFactory.newInstance();
			saxFact.setValidating(false);
			saxFact.setXIncludeAware(false);
			saxFact.setNamespaceAware(false);
			SAXParser sp = saxFact.newSAXParser();
			RNAMLParser RNAMLData = new RNAMLParser();
			sp.parse(new InputSource(r), RNAMLData);
			
			/*XMLReader xr = XMLReaderFactory.createXMLReader();
			RNAMLParser RNAMLData = new RNAMLParser();
			xr.setContentHandler(RNAMLData);
			xr.setErrorHandler(RNAMLData);
			xr.setEntityResolver(RNAMLData);
			xr.parse(new InputSource(r));*/

			setRNA(RNAMLData.getSequence(), RNAMLData.getBasicPlanarStructure());
			Vector<RNAMLParser.BPTemp> bps = RNAMLData.getPlanarBPs();
			for (int i = 0; i < bps.size(); i++) {
				RNAMLParser.BPTemp bp = bps.get(i);
				int bp5 = bp.pos5 - 1;
				int bp3 = bp.pos3 - 1;
				ModeleBase mb = _listeBases.get(bp5);
				ModeleBase part = _listeBases.get(bp3);
				ModeleStyleBP newStyle = bp.createBPStyle(mb, part);
				mb.setStyleBP(newStyle);
				part.setStyleBP(newStyle);
			}

			Vector<Integer> basenumbers = RNAMLData.getBaseNumbers();
			for (int i = 0; i < _listeBases.size(); i++) {
				ModeleBase n =  _listeBases.get(i);
				n.setBaseNumber(basenumbers.get(i));
			}

			Vector<RNAMLParser.BPTemp> bpsAux = RNAMLData.getAuxBPs();
			for (int i = 0; i < bpsAux.size(); i++) {
				RNAMLParser.BPTemp bp = bpsAux.get(i);
				int bp5 = bp.pos5 - 1;
				int bp3 = bp.pos3 - 1;
				
				//bp.createBPStyle(mb5, mb3);
				//addBPAux(bp5,bp3);
				ModeleBase mb = _listeBases.get(bp5);
				ModeleBase part = _listeBases.get(bp3);
				ModeleStyleBP newStyle = bp.createBPStyle(mb, part);
				addBPAux(bp5,bp3,newStyle);
			}
			loadOk = true;

		} catch (IOException ioe) {
			throw new ExceptionLoadingFailed(
					"Couldn't load file due to I/O or security policy issues.",
					"");
		} catch (Exception ge) {
          ge.printStackTrace();
		}
		return loadOk;
	}

	private boolean loadSecStrDBN(Reader r) throws ExceptionLoadingFailed,
			ExceptionPermissionDenied, ExceptionUnmatchedClosingParentheses,
			ExceptionFileFormatOrSyntax {
		boolean loadOk = false;
		try {
			BufferedReader fr = new BufferedReader(r);
			String line = fr.readLine();
			String title = "";
			String seqTmp = "";
			String strTmp = "";
			while ((line != null) && (strTmp.equals(""))) {
				line = line.trim();
				if (!line.startsWith(">")) {
					if (seqTmp.equals("")) {
						seqTmp = line;
					} else {
						strTmp = line;
					}
				}
				else
				{
					title = line.substring(1).trim();
				}
				line = fr.readLine();
			}
			if (strTmp.length() != 0) {
				setRNA(seqTmp, strTmp);
				this.setName(title);
				loadOk = true;
			}
		} catch (IOException e) {
			throw new ExceptionLoadingFailed(e.getMessage(), "");
		}
		return loadOk;
	}

	public void loadSecStr(Reader r) throws ExceptionFileFormatOrSyntax {
		loadSecStr(r,FILE_TYPE_UNKNOWN);
		}

	public void loadSecStr(Reader r, int fileType) throws ExceptionFileFormatOrSyntax 
	{			
		switch(fileType)
		{		
			case (FILE_TYPE_DBN):
			{
				try {
					boolean ok = loadSecStrDBN(r);
					if (ok) return;
				} catch (Exception e) { }
			}
			break;
			case (FILE_TYPE_CT):
			{
				try {
					boolean ok = loadSecStrCT(r);
					if (ok) return;
				} catch (Exception e) { }
			}
			break;
			case (FILE_TYPE_BPSEQ):
			{
				try {
					boolean ok = loadSecStrBPSEQ(r);
					if (ok) return;
				} catch (Exception e) { }
			}
			break;
			case (FILE_TYPE_RNAML):
			{
				try {
					boolean ok = loadSecStrRNAML(r);
					if (ok) return;
				} catch (Exception e) { }
			}
			break;
			case (FILE_TYPE_UNKNOWN):
			{
				BufferedReader buf = new BufferedReader(r);

				try {
					buf.mark(1000000);
					try {
						boolean ok = loadSecStrCT(buf);
						if (ok)
							return;
					} catch (Exception e) {
					}
					buf.reset();
					try {
						boolean ok = loadSecStrBPSEQ(buf);
						if (ok)
							return;
					} catch (Exception e) {
					}
					buf.reset();
					try {
						boolean ok = loadSecStrDBN(buf);
						if (ok)
							return;
					} catch (Exception e) {
						e.printStackTrace();
					}
					buf.reset();
					try {
						boolean ok = loadSecStrRNAML(buf);
						if (ok)
							return;
					} catch (ExceptionLoadingFailed e2)
					{
						e2.printStackTrace();
					} catch (Exception e) {
						e.printStackTrace();
					}
					try {
						boolean ok = loadSecStrRNAML(buf);
						if (ok)
							return;
					} catch (ExceptionLoadingFailed e2)
					{
						e2.printStackTrace();
					} catch (Exception e) {
						e.printStackTrace();
					}
					buf.reset();
				} catch (IOException e2) {
					e2.printStackTrace();
				}
			}
		}		
		throw new ExceptionFileFormatOrSyntax("");
	}

	public static int guessFileTypeFromExtension(String path)
	{
		if (path.toLowerCase().endsWith("ml"))
		{ return RNA.FILE_TYPE_RNAML; }
		else if (path.toLowerCase().endsWith("dbn")||path.toLowerCase().endsWith("faa"))
		{ return RNA.FILE_TYPE_DBN; }
		else if (path.toLowerCase().endsWith("ct"))
		{ return RNA.FILE_TYPE_CT; }
		else if (path.toLowerCase().endsWith("bpseq"))
		{ return RNA.FILE_TYPE_BPSEQ; }
		
		return RNA.FILE_TYPE_UNKNOWN; 			

	}
	
	public void loadSecStr(String path) throws ExceptionExportFailed,
			ExceptionPermissionDenied, ExceptionLoadingFailed,
			ExceptionFileFormatOrSyntax, ExceptionUnmatchedClosingParentheses,
			FileNotFoundException {
		FileReader fr = null;
		try {
			fr = new FileReader(path); 
			int type = guessFileTypeFromExtension(path);
			loadSecStr(fr,type);
		} catch (ExceptionFileFormatOrSyntax e) {
			if (fr != null)
				try {fr.close();} catch(IOException e2){}
			e.setPath(path);
			throw e;
		}
	}

	public void set_listeBases(ArrayList<ModeleBase> _liste) {
		this._listeBases = _liste;
	}

	public void addVARNAListener(InterfaceVARNAListener rl) {
		_listeVARNAListener.add(rl);
	}

	public void warningEmition(String warningMessage) {
		for (int i = 0; i < _listeVARNAListener.size(); i++) {
			_listeVARNAListener.get(i).warningEmitted(warningMessage);
		}
	}

	public void applyStyleOnBases(ArrayList<Integer> basesList,
			ModeleStyleBase style) {
		for (int i = 1; i < basesList.size(); i++) {
			_listeBases.get(basesList.get(i)).setStyleBase(style);
		}
	}

	private int[] correctReciprocity(int[] str) {
		int[] result = new int[str.length];
		for (int i = 0; i < str.length; i++) {
			if (str[i] != -1) {
				if (i == str[str[i]]) {
					result[i] = str[i];
				}
				else {
					str[str[i]] = i;
				}
			}
			else
			{
				result[i] = -1;
			}
		}
		return result;
	}

	private void applyStruct(int[] str) throws ExceptionFileFormatOrSyntax {
		str = correctReciprocity(str);

		int[] planarSubset = RNAMLParser.planarize(str);
		_structureAux.clear();

		for (int i = 0; i < planarSubset.length; i++) {
			_listeBases.get(i).setElementStructure(planarSubset[i]);
			if (planarSubset[i] > i) {
                addBP(i,planarSubset[i]);
			}
			else if ((planarSubset[i]!=str[i]))
			{
				addBPAux(i,str[i]);
			}
		}

	}

	@SuppressWarnings("unused")
	private void initStruct(int size) {
		for (int i = 0; i < size; i++) {
			_listeBases.get(i).setElementStructure(-1);
		}
	}

	public ArrayList<ModeleBase> get_listeBases() {
		return _listeBases;
	}

	public boolean is_comparisonMode() {
		return _comparisonMode;
	}

	public void set_comparisonMode(boolean mode) {
		_comparisonMode = mode;
	}

	public int getSize() {
		return _listeBases.size();
	}

	public ArrayList<Integer> findAll() {
		ArrayList<Integer> listAll = new ArrayList<Integer>();
		for (int i = 0; i < get_listeBases().size(); i++) {
			listAll.add(i);
		}
		return listAll;
	}

	public ArrayList<Integer> findBulge(int index) {
		ArrayList<Integer> listUp = new ArrayList<Integer>();
		if (get_listeBases().get(index).getElementStructure() == -1) {
			int i = index;
			boolean over = false;
			while ((i < get_listeBases().size()) && !over) {
				int j = get_listeBases().get(i).getElementStructure();
				if (j == -1) {
					listUp.add(i);
					i++;
				} else {
					over = true;
				}
			}
			i = index - 1;
			over = false;
			while ((i >= 0) && !over) {
				int j = get_listeBases().get(i).getElementStructure();
				if (j == -1) {
					listUp.add(i);
					i--;
				} else {
					over = true;
				}
			}
		}
		return listUp;
	}

	public ArrayList<Integer> findStem(int index) {
		ArrayList<Integer> listUp = new ArrayList<Integer>();
		int i = index;
		do {
			listUp.add(i);
			int j = get_listeBases().get(i).getElementStructure();
			if (j == -1) {
				i = (i + 1) % getSize();
			} else {
				if ((j < i) && (index <= i) && (j <= index)) {
					i = j;
				} else {
					i = (i + 1) % getSize();
				}
			}
		} while (i != index);
		return listUp;
	}

	public int getHelixCountOnLoop(int indice) {
		int cptHelice = 0;
		if (indice < 0 || indice >= get_listeBases().size())
			return cptHelice;
		int i = indice;
		int j = get_listeBases().get(i).getElementStructure();
		// Only way to distinguish "supporting base-pair" from others
		boolean justJumped = false;
		if ((j != -1) && (j < i)) {
			i = j + 1;
			indice = i;
		}
		do {
			j = get_listeBases().get(i).getElementStructure();
			if ((j != -1) && (!justJumped)) {
				i = j;
				justJumped = true;
				cptHelice++;
			} else {
				i = (i + 1) % get_listeBases().size();
				justJumped = false;
			}
		} while (i != indice);
		return cptHelice;
	}

	public ArrayList<Integer> findLoop(int indice) {
		return findLoopForward(indice);
	}

	public ArrayList<Integer> findLoopForward(int indice) {
		ArrayList<Integer> base = new ArrayList<Integer>();
		if (indice < 0 || indice >= get_listeBases().size())
			return base;
		int i = indice;
		int j = get_listeBases().get(i).getElementStructure();
		// Only way to distinguish "supporting base-pair" from others
		boolean justJumped = false;
		if (j != -1) {
			i = Math.min(i, j) + 1;
			indice = i;
		}
		do {
			base.add(i);
			j = get_listeBases().get(i).getElementStructure();
			if ((j != -1) && (!justJumped)) {
				i = j;
				justJumped = true;
			} else {
				i = (i + 1) % get_listeBases().size();
				justJumped = false;
			}
		} while (i != indice);
		return base;
	}

	public ArrayList<Integer> findPair(int indice) {
		ArrayList<Integer> base = new ArrayList<Integer>();
		int j = get_listeBases().get(indice).getElementStructure();
		if (j != -1) {
			base.add(Math.min(indice, j));
			base.add(Math.max(indice, j));
		}

		return base;

	}

	public ArrayList<Integer> findLoopBackward(int indice) {
		ArrayList<Integer> base = new ArrayList<Integer>();
		if (indice < 0 || indice >= get_listeBases().size())
			return base;
		int i = indice;
		int j = get_listeBases().get(i).getElementStructure();
		// Only way to distinguish "supporting base-pair" from others
		boolean justJumped = false;
		if (j != -1) {
			i = Math.min(i, j) - 1;
			indice = i;
		}
		if (i < 0) {
			return base;
		}
		do {
			base.add(i);
			j = get_listeBases().get(i).getElementStructure();
			if ((j != -1) && (!justJumped)) {
				i = j;
				justJumped = true;
			} else {
				i = (i + get_listeBases().size() - 1) % get_listeBases().size();
				justJumped = false;
			}
		} while (i != indice);
		return base;
	}

	public ArrayList<Integer> findHelix(int indice) {
		ArrayList<Integer> list = new ArrayList<Integer>();
		if (get_listeBases().get(indice).getElementStructure() != -1) {
			list.add(indice);
			list.add(get_listeBases().get(indice).getElementStructure());
			int i = 1, prec = get_listeBases().get(indice)
					.getElementStructure();
			while (indice + i < get_listeBases().size()
					&& get_listeBases().get(indice + i).getElementStructure() != -1
					&& get_listeBases().get(indice + i).getElementStructure() == prec - 1) {
				list.add(indice + i);
				list.add(get_listeBases().get(indice + i)
						.getElementStructure());
				prec = get_listeBases().get(indice + i).getElementStructure();
				i++;
			}
			i = -1;
			prec = get_listeBases().get(indice).getElementStructure();
			while (indice + i >= 0
					&& get_listeBases().get(indice + i).getElementStructure() != -1
					&& get_listeBases().get(indice + i).getElementStructure() == prec + 1) {
				list.add(indice + i);
				list.add(get_listeBases().get(indice + i)
						.getElementStructure());
				prec = get_listeBases().get(indice + i).getElementStructure();
				i--;
			}
		}
		return list;
	}

	public ArrayList<Integer> find3Prime(int indice) {
		ArrayList<Integer> list = new ArrayList<Integer>();
		boolean over = false;
		while ((indice >= 0) && !over) {
			over = (get_listeBases().get(indice).getElementStructure() != -1);
			indice--;
		}
		indice++;
		if (over) {
			indice++;
		}
		for (int i = indice; i < get_listeBases().size(); i++) {
			list.add(i);
			if (get_listeBases().get(i).getElementStructure() != -1) {
				return new ArrayList<Integer>();
			}
		}
		return list;
	}

	public ArrayList<Integer> find5Prime(int indice) {
		ArrayList<Integer> list = new ArrayList<Integer>();
		for (int i = 0; i <= indice; i++) {
			list.add(i);
			if (get_listeBases().get(i).getElementStructure() != -1) {
				return new ArrayList<Integer>();
			}
		}
		return list;
	}

	public Double get_spaceBetweenBases() {
		return _spaceBetweenBases;
	}

	public void set_spaceBetweenBases(Double betweenBases) {
		_spaceBetweenBases = betweenBases;
	}

	public static Double angle(Point2D.Double p1, Point2D.Double p2,
			Point2D.Double p3) {
		Double alpha = Math.atan2(p1.y - p2.y, p1.x - p2.x);
		Double beta = Math.atan2(p3.y - p2.y, p3.x - p2.x);
		Double angle = (beta - alpha);

		// Correction de l'angle pour le resituer entre 0 et 2PI
		while (angle < 0.0 || angle > 2 * Math.PI) {
			if (angle < 0.0)
				angle += 2 * Math.PI;
			else if (angle > 2 * Math.PI)
				angle -= 2 * Math.PI;
		}
		return angle;
	}

	public ArrayList<Integer> findNonPairedBaseGroup(Integer get_nearestBase) {
		// detection 3', 5', bulge
		ArrayList<Integer> list = new ArrayList<Integer>();
		int indice = get_nearestBase;
		boolean nonpairedUp = true, nonpairedDown = true;
		while (indice < get_listeBases().size() && nonpairedUp) {
			if (get_listeBases().get(indice).getElementStructure() == -1) {
				list.add(indice);
				indice++;
			} else {
				nonpairedUp = false;
			}
		}
		indice = get_nearestBase - 1;
		while (indice >= 0 && nonpairedDown) {
			if (get_listeBases().get(indice).getElementStructure() == -1) {
				list.add(indice);
				indice--;
			} else {
				nonpairedDown = false;
			}
		}
		return list;
	}

	public boolean getDrawn() {
		return _drawn;
	}

	public ArrayList<ModeleStyleBP> getStructureAux() {
		return _structureAux;
	}

	public int getIndexFromBaseNumber(int num) {
		for (int i = 0; i < this._listeBases.size(); i++) {
			if (_listeBases.get(i).getBaseNumber() == num) {
				return i;
			}
		}
		return -1;
	}

	/**
	 * Adds a base pair to this RNA's structure. Tries to add it to the
	 * secondary structure first, eventually adding it to the 'tertiary'
	 * interactions if it clashes with the current secondary structure.
	 * 
	 * @param a
	 *            - Base number of the origin of this base pair
	 * @param b
	 *            - Base number of the destination of this base pair
	 */

	public void addBPToStructure(int a, int b) {
		int i = getIndexFromBaseNumber(a);
		int j = getIndexFromBaseNumber(b);
		ModeleBase part5 = _listeBases.get(i);
		ModeleBase part3 = _listeBases.get(j);
		ModeleStyleBP msbp = new ModeleStyleBP(part5, part3);
		addBPToStructure(a, b, msbp);
	}

	/**
	 * Adds a base pair to this RNA's structure. Tries to add it to the
	 * secondary structure first, possibly adding it to the 'tertiary'
	 * interactions if it clashes with the current secondary structure.
	 * 
	 * @param a
	 *            - Base number of the origin of this base pair
	 * @param b
	 *            - Base number of the destination of this base pair
	 */

	public void addBPToStructure(int a, int b, ModeleStyleBP msbp) {
		int i = getIndexFromBaseNumber(a);
		int j = getIndexFromBaseNumber(b);
		if (j < i) {
			int k = j;
			j = i;
			i = k;
		}
		if (i != -1) {
			if ((_listeBases.get(i).getElementStructure() != -1)
					|| (_listeBases.get(j).getElementStructure() != -1)) {
				addBPAux(i, j, msbp);
				return;
			}

			for (int k = i + 1; k < j; k++) {
				ModeleBase tmp = _listeBases.get(k);
				int l = tmp.getElementStructure();
				if (l != -1) {
					if ((l <= i) || (l >= j)) {
						addBPAux(i, j, msbp);
						return;
					}
				}
			}
			addBP(i, j, msbp);
		}
	}

	public void addBP(int i, int j) {
		if (j < i) {
			int k = j;
			j = i;
			i = k;
		}

		ModeleBase part5 = _listeBases.get(i);
		ModeleBase part3 = _listeBases.get(j);
		ModeleStyleBP msbp = new ModeleStyleBP(part5, part3);
		addBP(i, j, msbp);
	}

	public void addBP(int i, int j, ModeleStyleBP msbp) {
		if (j < i) {
			int k = j;
			j = i;
			i = k;
		}
		ModeleBase part5 = _listeBases.get(i);
		ModeleBase part3 = _listeBases.get(j);
		msbp.setPartner5(part5);
		msbp.setPartner3(part3);
		part5.setElementStructure(j, msbp);
		part3.setElementStructure(i, msbp);
	}

	public void addBPAux(int i, int j) {
		ModeleBase part5 = _listeBases.get(i);
		ModeleBase part3 = _listeBases.get(j);
		ModeleStyleBP msbp = new ModeleStyleBP(part5, part3);
		addBPAux(i, j, msbp);
	}

	public void addBPAux(int i, int j, ModeleStyleBP msbp) {
		if (j < i) {
			int k = j;
			j = i;
			i = k;
		}
		ModeleBase part5 = _listeBases.get(i);
		ModeleBase part3 = _listeBases.get(j);
		msbp.setPartner5(part5);
		msbp.setPartner3(part3);
		_structureAux.add(msbp);
	}

	public ModeleStyleBP getBPStyle(int i, int j) {
		ModeleStyleBP result = null;
		if (i > j) {
			int k = j;
			j = i;
			i = k;
		}
		if (_listeBases.get(i).getElementStructure() == j) {
			result = _listeBases.get(i).getStyleBP();
		}
		for (int k = 0; k < _structureAux.size(); k++) {
			ModeleStyleBP bp = _structureAux.get(k);
			if ((bp.getPartner5().getIndex() == i)
					&& (bp.getPartner3().getIndex() == j)) {
				result = bp;
			}
		}
		return result;
	}

	public ArrayList<ModeleStyleBP> getAuxBPs(int i)
	{
		ArrayList<ModeleStyleBP> result = new ArrayList<ModeleStyleBP>(); 
		for (ModeleStyleBP bp : _structureAux) 
		{
			if ((bp.getPartner5().getIndex() == i) || (bp.getPartner3().getIndex() == i)) {
				result.add(bp);
			}
		}
		return result;
	}
	
	public void setBaseInnerColor(Color c) {
		for (int i = 0; i < _listeBases.size(); i++) {
			ModeleBase mb = _listeBases.get(i);
			mb.getStyleBase().set_base_inner_color(c);
		}
	}

	public void setBaseNumbersColor(Color c) {
		for (int i = 0; i < _listeBases.size(); i++) {
			ModeleBase mb = _listeBases.get(i);
			mb.getStyleBase().set_base_number_color(c);
		}
	}

	public void setBaseNameColor(Color c) {
		for (int i = 0; i < _listeBases.size(); i++) {
			ModeleBase mb = _listeBases.get(i);
			mb.getStyleBase().set_base_name_color(c);
		}
	}

	public void setBaseOutlineColor(Color c) {
		for (int i = 0; i < _listeBases.size(); i++) {
			ModeleBase mb = _listeBases.get(i);
			mb.getStyleBase().set_base_outline_color(c);
		}
	}

	public boolean getFlatExteriorLoop()
	{
		return _flatExteriorLoop;
	}

	public void setFlatExteriorLoop(boolean on)
	{
		_flatExteriorLoop = on;
	}
	
	public String getName()
	{
		return _name;
	}

	public void setName(String n)
	{
		_name = n;
	}
	
	public ArrayList<TextAnnotation> getAnnotations()
	{
		return _listeAnnotations;
	}

	public boolean removeAnnotation(TextAnnotation t)
	{
		return _listeAnnotations.remove(t);
	}

	public void addAnnotation(TextAnnotation t)
	{
		_listeAnnotations.add(t);
	}

	public void removeAnnotation(String filter)
	{
		ArrayList<TextAnnotation> condamne = new ArrayList<TextAnnotation>();
		for (TextAnnotation t : _listeAnnotations)
		{
			if (t.getTexte().contains(filter))
			{
				condamne.add(t);				
			}
		}
		for (TextAnnotation t : condamne)
		{
			_listeAnnotations.remove(t);
		}
	}

	public void clearAnnotations()
	{
		_listeAnnotations.clear();
	}
	
	private boolean _strandEndsAnnotated = false;
	
	public void autoAnnotateStrandEnds()
	{
		if (! _strandEndsAnnotated)
		{
		int tailleListBases=_listeBases.size();
		boolean endAnnotate =false;
		addAnnotation(new TextAnnotation("5'", _listeBases.get(0)));
		for(int i=0; i<_listeBases.size()-1; i++){
			int realposA = _listeBases.get(i).getBaseNumber();
			int realposB = _listeBases.get(i+1).getBaseNumber();
			if(realposB-realposA!=1){
				addAnnotation(new TextAnnotation("3'", _listeBases.get(i)));
				addAnnotation(new TextAnnotation("5'", _listeBases.get(i+1)));
				if(i+1==_listeBases.size()-1){
					endAnnotate=true;
				}
			}
		}
		if(!endAnnotate){
			addAnnotation(new TextAnnotation("3'", _listeBases.get(tailleListBases-1)));
		}
		_strandEndsAnnotated = true;
		}
		else
		{
			removeAnnotation("3'");
			removeAnnotation("5'");
			_strandEndsAnnotated = false;
		}
	}
	
	public void autoAnnotateHelices()
	{
		Stack<Integer> p = new Stack<Integer>();
		p.push(0);
		int nbH = 1;
		while(!p.empty())
		{
			int i = p.pop();
			if (i<_listeBases.size())
			{
				ModeleBase mb = _listeBases.get(i); 
				int j = mb.getElementStructure(); 
				if (j==-1)
				{ p.push(i+1);	}
				else
				{
					if (j>i)
					{
						ModeleBase mbp = _listeBases.get(j); 
						p.push(j+1);
						ArrayList<ModeleBase> h = new ArrayList<ModeleBase>(); 
						int k = 1;
						while(mb.getElementStructure()==mbp.getIndex())
						{
							h.add(mb);
							h.add(mbp);
							mb = _listeBases.get(i+k);
							mbp = _listeBases.get(j-k); 
						
							k++;
						}
						try {
							addAnnotation(new TextAnnotation("H"+nbH++, h , TextAnnotation.HELIX));
						} catch (Exception e) {
							e.printStackTrace();
						}
						p.push(i+k);
					}
				}
			}
		}
	}

	public void autoAnnotateTerminalLoops()
	{
		Stack<Integer> p = new Stack<Integer>();
		p.push(0);
		int nbT = 1;
		while(!p.empty())
		{
			int i = p.pop();
			if (i<_listeBases.size())
			{
			ModeleBase mb = _listeBases.get(i); 
			int j = mb.getElementStructure(); 
			if (j==-1)
			{ 	
				int k = 1;
				ArrayList<ModeleBase> t = new ArrayList<ModeleBase>();
				while ((i+k<getSize())&&(mb.getElementStructure()==-1))
				{
					t.add(mb);
					mb = _listeBases.get(i+k);
					k++;
				}
				if (mb.getElementStructure()!=-1)
				{
					if (mb.getElementStructure()==i-1)
					{
						try {
							t.add(_listeBases.get(i-1));
							t.add(_listeBases.get(i+k-1));
							addAnnotation(new TextAnnotation("T"+nbT++, t , TextAnnotation.LOOP));
						} catch (Exception e) {
							e.printStackTrace();
						}
					}
					p.push(i+k-1);
				}
				
			}
			else
			{
				if (j>i)
				{
					p.push(j+1);
					p.push(i+1);
				}
			}
			}
		}
	}

	public void autoAnnotateInteriorLoops()
	{
		Stack<Integer> p = new Stack<Integer>();
		p.push(0);
		int nbT = 1;
		while(!p.empty())
		{
			int i = p.pop();
			if (i<_listeBases.size())
			{
				ModeleBase mb = _listeBases.get(i); 
				int j = mb.getElementStructure(); 
				if (j==-1)
				{ 	
					int k = i+1;
					ArrayList<ModeleBase> t = new ArrayList<ModeleBase>();
					boolean terminal = true;
					while ((k<getSize())&&((mb.getElementStructure()>=i)||(mb.getElementStructure()==-1)))
					{
						t.add(mb);
						mb = _listeBases.get(k);
						if ((mb.getElementStructure()==-1)||(mb.getElementStructure()<k))
							k++;
						else 
						{
							p.push(k);
							terminal = false;
							k = mb.getElementStructure();
						}
					}
					if (mb.getElementStructure()!=-1)
					{
						if ((mb.getElementStructure()==i-1)&& !terminal)
						{
							try {
								t.add(_listeBases.get(i-1));
								t.add(_listeBases.get(k-1));
								addAnnotation(new TextAnnotation("I"+nbT++, t , TextAnnotation.LOOP));
							} catch (Exception e) {
								e.printStackTrace();
							}
							p.push(k-1);
						}
					}	
				}
				else
				{
				if (j>i)
				{
					p.push(i+1);
				}
			}
		}
		}
	}
	
	@SuppressWarnings("unchecked")
	public TextAnnotation getAnnotation(int type, ModeleBase base)
	{
		TextAnnotation result = null;
		for(TextAnnotation t : _listeAnnotations)
		{
			if (t.getType()==type)
			{
				switch(type)
				{
				    case(TextAnnotation.BASE):
				    	if (base == (ModeleBase) t.getAncrage())
				    	  return t;
					break;
				    case(TextAnnotation.HELIX):
				    case(TextAnnotation.LOOP):
				    {
				    	ArrayList<ModeleBase> mbl = (ArrayList<ModeleBase>) t.getAncrage();
				    	if (mbl.contains(base))
				    	  return t;
				    }
					break;
				}
			}
		}
		return result;
	}
	
	private ArrayList<ChemProbAnnotation> _ChemProbAnnotations = new ArrayList<ChemProbAnnotation>();
	
	public void addChemProbAnnotation(ChemProbAnnotation cpa)
	{
		_ChemProbAnnotations.add(cpa);
	}
	
	public ArrayList<ChemProbAnnotation> getChemProbAnnotations()
	{
		return _ChemProbAnnotations;
	}

	public void setColorMapValues(Double[] values, ModeleColorMap cm)
	{
		setColorMapValues(values, cm, false);
	}
	
	public void adaptColorMapToValues(ModeleColorMap cm)
	{
		double min = Double.MAX_VALUE;
		double max = Double.MIN_VALUE;
		for (int i=0;i<Math.min(_listeBases.size(),_listeBases.size());i++)
		{
			ModeleBase mb = _listeBases.get(i);
			max = Math.max(max,mb.getValue());
			min = Math.min(min,mb.getValue());
		}
		cm.rescale(min, max);		
	}
	
	  public void readValues(Reader r, ModeleColorMap cm)
	  {
		  try {
			  StreamTokenizer st = new StreamTokenizer(r);
			  st.eolIsSignificant(true);
			  ArrayList<ArrayList<Double>> vals = new ArrayList<ArrayList<Double>>();
			  ArrayList<Double> curVals = new ArrayList<Double>(); 
			  int type = st.nextToken();
			  while (type != StreamTokenizer.TT_EOF)
			  {
				  switch(type)
				  {
				  case (StreamTokenizer.TT_NUMBER): 
					  curVals.add(st.nval);
				  break;
				  case (StreamTokenizer.TT_EOL):
					  if (curVals.size()>0)
					  {
						  vals.add(curVals);
						  curVals = new ArrayList<Double>();
					  }
				  break;
				  }
				  type = st.nextToken();
			  }
			  if (curVals.size()>0) 
				  vals.add(curVals);
			  
			  Double[] v = new Double[vals.size()];
			  for (int i=0;i<Math.min(vals.size(),getSize());i++)
			  {
				  ArrayList<Double> tab = vals.get(i);
				  v[i] = tab.get(tab.size()-1);
			  }
			  setColorMapValues(v, cm,  true);
		  } catch (IOException e) {
			  e.printStackTrace();
		  }
	  }


	public void setColorMapValues(Double[] values, ModeleColorMap cm, boolean rescaleColorMap)
	{
		if (values.length>0)
		{
			for (int i=0;i<Math.min(values.length,_listeBases.size());i++)
			{
				ModeleBase mb = _listeBases.get(i);
				mb.setValue(values[i]);
			}
			if (rescaleColorMap){
				adaptColorMapToValues(cm);
			}
		}
	}

	public Double[] getColorMapValues()
	{
		Double[] values = new Double[_listeBases.size()];
		for (int i=0;i<_listeBases.size();i++)
		{
			values[i] = _listeBases.get(i).getValue();
		}
		return values;
	}
	
	public void rescaleColorMap(ModeleColorMap cm)
	{
		Double max = Double.MIN_VALUE;
		Double min = Double.MAX_VALUE;
		for (int i=0;i<_listeBases.size();i++)
		{
			Double value = _listeBases.get(i).getValue();
			max = Math.max(max, value);
			min = Math.min(min, value);
		}
		cm.rescale(min, max);
	}
	
	public void setColorMapValue(int index, double value, ModeleColorMap cm)
	{
		Double[] values = new Double[1];
		values[0] = value;
		setColorMapValues(values, cm, false);
	}
	
	public void setSequence(String s)
	{
		int i = 0;
		int j = 0;
		while ((i<s.length())&&(j<_listeBases.size()))
		{
			ModeleBase mb = _listeBases.get(j);
			if (mb instanceof ModeleBaseNucleotide)
			{
				((ModeleBaseNucleotide)mb).set_c(s.charAt(i));
				i++; j++;
			}
			else if (mb instanceof ModeleBasesComparison)
			{
				((ModeleBasesComparison)mb).set_base1(s.charAt(i));
				((ModeleBasesComparison)mb).set_base2(s.charAt(i+1));
				i+=2;j++;
			}
			else j++;
 		}
	}
	
	public void eraseSequence()
	{
		int j = 0;
		while ((j<_listeBases.size()))
		{
			ModeleBase mb = _listeBases.get(j);
			if (mb instanceof ModeleBaseNucleotide)
			{
				((ModeleBaseNucleotide)mb).set_c(' ');
				j++;
			}
			else if (mb instanceof ModeleBasesComparison)
			{
				((ModeleBasesComparison)mb).set_base1(' ');
				((ModeleBasesComparison)mb).set_base2(' ');
				j++;
			}
			else j++;
 		}		
	}
	
	
    public RNA clone ()
    {
        try
        {
            ByteArrayOutputStream out = new ByteArrayOutputStream ();
            ObjectOutputStream oout = new ObjectOutputStream (out);
            oout.writeObject (this);
            
            ObjectInputStream in = new ObjectInputStream (
                new ByteArrayInputStream (out.toByteArray ()));
            return (RNA)in.readObject ();
        }
        catch (Exception e)
        {
            throw new RuntimeException ("cannot clone class [" +
                this.getClass ().getName () + "] via serialization: " +
                e.toString ());
        }
    }
    
    public ModeleBase getBaseAt(int index)
    {
    	return this._listeBases.get(index);
    }

    public ArrayList<ModeleBase> getBasesAt(Collection<? extends Integer> indices)
    {
    	ArrayList<ModeleBase> mbs = new ArrayList<ModeleBase>();
    	for (int i: indices)
    	{
    		mbs.add(getBaseAt(i));
    	}
    	return mbs;
    }

    public ArrayList<ModeleBase> getBasesBetween(int from, int to)
    {
    	ArrayList<ModeleBase> mbs = new ArrayList<ModeleBase>();
    	int bck = Math.min(from, to);
    	to = Math.max(from, to);
    	from = bck;
    	for (int i=from;i<=to;i++)
    	{  mbs.add(getBaseAt(i)); }
    	return mbs;
    }

    public void addHighlightRegion(HighlightRegionAnnotation n)
    {
    	_listeRegionHighlights.add(n);    	
    }

	public void removeHighlightRegion(HighlightRegionAnnotation n)
	{
		_listeRegionHighlights.remove(n);
	}
   
	public void removeChemProbAnnotation(ChemProbAnnotation a)
	{
		_ChemProbAnnotations.remove(a);
	}
	
    public void addHighlightRegion(int from, int to, Color fill, Color outline, double radius)
    {
    	_listeRegionHighlights.add(new HighlightRegionAnnotation(getBasesBetween(from, to),fill,outline,radius));    	
    }
    
    
    public void addHighlightRegion(int from, int to)
    {
    	_listeRegionHighlights.add(new HighlightRegionAnnotation(getBasesBetween(from, to)));    	
    }

    public ArrayList<HighlightRegionAnnotation> getHighlightRegion()
    {
    	return _listeRegionHighlights;    	
    }


	/**
	 * Rotates the RNA coordinates by a certain angle
	 * 
	 * @param angleDegres
	 *            Rotation angle, in degrees
	 */
	public void globalRotation(Double angleDegres) {
		if (_listeBases.size() > 0) {

			// angle en radian
			Double angle = angleDegres * Math.PI / 180;

			// initialisation du minimum et dumaximum
			Double maxX = _listeBases.get(0).getCoords().x;
			Double maxY = _listeBases.get(0).getCoords().y;
			Double minX = _listeBases.get(0).getCoords().x;
			Double minY = _listeBases.get(0).getCoords().y;
			// mise a jour du minimum et du maximum
			for (int i = 0; i < _listeBases.size(); i++) {
				if (_listeBases.get(i).getCoords().getX() < minX)
					minX = _listeBases.get(i).getCoords().getX();
				if (_listeBases.get(i).getCoords().getY() < minY)
					minY = _listeBases.get(i).getCoords().getY();
				if (_listeBases.get(i).getCoords().getX() > maxX)
					maxX = _listeBases.get(i).getCoords().getX();
				if (_listeBases.get(i).getCoords().getX() > maxY)
					maxY = _listeBases.get(i).getCoords().getY();
			}
			// creation du point central
			Point2D.Double centre = new Point2D.Double((maxX - minX) / 2,
					(maxY - minY) / 2);
			Double x, y;
			for (int i = 0; i < _listeBases.size(); i++) {
				// application de la rotation au centre de chaque base
				// x' = cos(theta)*(x-xc) - sin(theta)*(y-yc) + xc
				x = Math.cos(angle)
						* (_listeBases.get(i).getCenter().getX() - centre.x)
						- Math.sin(angle)
						* (_listeBases.get(i).getCenter().getY() - centre.y)
						+ centre.x;
				// y' = sin(theta)*(x-xc) + cos(theta)*(y-yc) + yc
				y = Math.sin(angle)
						* (_listeBases.get(i).getCenter().getX() - centre.x)
						+ Math.cos(angle)
						* (_listeBases.get(i).getCenter().getY() - centre.y)
						+ centre.y;
				_listeBases.get(i).setCenter(
						new Point2D.Double(x, y));

				// application de la rotation au coordonnees de chaque
				// base
				// x' = cos(theta)*(x-xc) - sin(theta)*(y-yc) + xc
				x = Math.cos(angle)
						* (_listeBases.get(i).getCoords().getX() - centre.x)
						- Math.sin(angle)
						* (_listeBases.get(i).getCoords().getY() - centre.y)
						+ centre.x;
				// y' = sin(theta)*(x-xc) + cos(theta)*(y-yc) + yc
				y = Math.sin(angle)
						* (_listeBases.get(i).getCoords().getX() - centre.x)
						+ Math.cos(angle)
						* (_listeBases.get(i).getCoords().getY() - centre.y)
						+ centre.y;
				_listeBases.get(i).setCoords(
						new Point2D.Double(x, y));
			}
		}
	}

}


