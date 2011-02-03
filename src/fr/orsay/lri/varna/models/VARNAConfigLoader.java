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
package fr.orsay.lri.varna.models;

/*
 * VARNA is a Java library for quick automated drawings RNA secondary structure
 * Copyright (C) 2007 Yann Ponty
 * 
 * This program is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 * 
 * You should have received a copy of the GNU General Public License along with
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.geom.Point2D;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URL;
import java.net.URLConnection;
import java.util.ArrayList;

import fr.orsay.lri.varna.VARNAPanel;
import fr.orsay.lri.varna.exceptions.ExceptionDrawingAlgorithm;
import fr.orsay.lri.varna.exceptions.ExceptionFileFormatOrSyntax;
import fr.orsay.lri.varna.exceptions.ExceptionLoadingFailed;
import fr.orsay.lri.varna.exceptions.ExceptionModeleStyleBaseSyntaxError;
import fr.orsay.lri.varna.exceptions.ExceptionNAViewAlgorithm;
import fr.orsay.lri.varna.exceptions.ExceptionNonEqualLength;
import fr.orsay.lri.varna.exceptions.ExceptionParameterError;
import fr.orsay.lri.varna.interfaces.InterfaceParameterLoader;
import fr.orsay.lri.varna.models.annotations.ChemProbAnnotation;
import fr.orsay.lri.varna.models.annotations.HighlightRegionAnnotation;
import fr.orsay.lri.varna.models.annotations.TextAnnotation;
import fr.orsay.lri.varna.models.rna.ModeleBase;
import fr.orsay.lri.varna.models.rna.ModeleColorMap;
import fr.orsay.lri.varna.models.rna.ModeleStyleBP;
import fr.orsay.lri.varna.models.rna.ModeleStyleBase;
import fr.orsay.lri.varna.models.rna.RNA;

/**
 * An RNA 2d Panel demo applet
 * 
 * @author Yann Ponty
 * 
 */

public class VARNAConfigLoader {

	private static final int MAXSTYLE = 50;

	// Applet Options
	public static String numColumnsOpt = "rows";
	public static String numRowsOpt = "columns";

	public static String algoOpt = "algorithm";
	public static String annotationsOpt = "annotations";
	public static String applyBasesStyleOpt = "applyBasesStyle";
	public static String auxBPsOpt = "auxBPs";
	public static String autoHelicesOpt = "autoHelices";
	public static String autoInteriorLoopsOpt = "autoInteriorLoops";
	public static String autoTerminalLoopsOpt = "autoTerminalLoops";

	public static String backboneColorOpt = "backbone";
	public static String backgroundColorOpt = "background";
	public static String baseInnerColorOpt = "baseInner";
	public static String baseNameColorOpt = "baseName";
	public static String baseNumbersColorOpt = "baseNum";
	public static String baseOutlineColorOpt = "baseOutline";
	public static String basesStyleOpt = "basesStyle";
	public static String borderOpt = "border";
	public static String bondColorOpt = "bp";
	public static String bpIncrementOpt = "bpIncrement";
	public static String bpStyleOpt = "bpStyle";

	public static String colorMapOpt = "colorMap";
	public static String colorMapCaptionOpt = "colorMapCaption";
	public static String colorMapDefOpt = "colorMapStyle";

	public static String chemProbOpt = "chemProb";
	public static String customBasesOpt = "customBases";
	public static String customBPsOpt = "customBPs";

	public static String drawNCOpt = "drawNC";
	public static String drawTertiaryOpt = "drawTertiary";
	public static String drawColorMapOpt = "drawColorMap";
	public static String drawBackboneOpt = "drawBackbone";

	public static String errorOpt = "error";

	public static String flatExteriorLoopOpt = "flat";

	public static String highlightRegionOpt = "highlightRegion";
	
	public static String nonStandardColorOpt = "nsBasesColor";

	public static String modifiableOpt = "modifiable";

	
	public static String periodNumOpt = "periodNum";

	public static String rotationOpt = "rotation";

	public static String sequenceOpt = "sequenceDBN";
	public static String structureOpt = "structureDBN";

	public static String titleOpt = "title";
	public static String titleColorOpt = "titleColor";
	public static String titleSizeOpt = "titleSize";

	public static String URLOpt = "url";

	public static String zoomOpt = "zoom";
	public static String zoomAmountOpt = "zoomAmount";

	public static String warningOpt = "warning";

	public static String comparisonModeOpt = "comparisonMode";
	public static String gapsBaseColorOpt = "gapsColor";
	public static String firstSequenceForComparisonOpt = "firstSequence";
	public static String secondSequenceForComparisonOpt = "secondSequence";
	public static String firstStructureForComparisonOpt = "firstStructure";
	public static String secondStructureForComparisonOpt = "secondStructure";

	// Applet assignable parameters
	private String _algo;
	public String _annotations;
	public String _chemProbs;
	private double _rotation;

	private String _sseq;
	private String _sstruct;

	private int _numRows;
	private int _numColumns;

	private String _title;
	private int _titleSize;
	private Color _titleColor;

	private String _auxBPs;
	private String _highlightRegion;
	
	private boolean _autoHelices;
	private boolean _autoInteriorLoops;
	private boolean _autoTerminalLoops;

	private boolean _drawBackbone;
	private Color _backboneColor;
	private Color _bondColor;
	private VARNAConfig.BP_STYLE _bpStyle;
	private Color _baseOutlineColor;
	private Color _baseInnerColor;
	private Color _baseNumColor;
	private Color _baseNameColor;
	private Color _gapsColor;
	private Color _nonStandardColor;

	private boolean _flatExteriorLoop;

	private String _customBases;
	private String _customBPs;

	private String _colorMapStyle;
	private String _colorMapCaption;
	private String _colorMapValues;

	private boolean _drawNC;
	private boolean _drawTertiary;
	private boolean _drawColorMap;

	private int _periodResNum;
	private Dimension _border;

	private Color _backgroundColor;

	private boolean _warning, _error;

	private boolean _modifiable;

	private double _zoom, _zoomAmount;

	private ArrayList<ModeleStyleBase> _basesStyleList;

	private boolean _comparisonMode;

	private String _firstSequence;
	private String _secondSequence;
	private String _firstStructure;
	private String _secondStructure;

	private VARNAPanel _mainSurface;

	private boolean _useNonStandardColor;
	private boolean _useGapsColor;
	private double _bpIncrement;

	private boolean _useInnerBaseColor;
	private boolean _useBaseNameColor;
	private boolean _useBaseNumbersColor;
	private boolean _useBaseOutlineColor;

	private String _URL;

	protected ArrayList<VARNAPanel> _VARNAPanelList = new ArrayList<VARNAPanel>();

	InterfaceParameterLoader _optionProducer;

	public VARNAConfigLoader(InterfaceParameterLoader il) {
		_optionProducer = il;
	}

	public ArrayList<VARNAPanel> createVARNAPanels()
			throws ExceptionParameterError,
			ExceptionModeleStyleBaseSyntaxError, ExceptionNonEqualLength,
			IOException, ExceptionFileFormatOrSyntax, ExceptionLoadingFailed {
		retrieveParametersValues();
		return _VARNAPanelList;
	}

	public int getNbRows() {
		return this._numRows;
	}

	public int getNbColumns() {
		return this._numColumns;
	}

	private void initValues() {

		// Applet assignable parameters
		_algo = "radiate";
		_auxBPs = "";
		_autoHelices = false;
		_autoInteriorLoops = false;
		_autoTerminalLoops = false;
		_annotations = "";
		_backgroundColor = VARNAConfig.DEFAULT_BACKGROUND_COLOR;
		_customBases = "";
		_customBPs = "";
		_chemProbs = "";
		
		
		_colorMapStyle = "";
		_colorMapValues = "";
		_colorMapCaption = "";
		_drawColorMap = false;
		
		_drawNC = true;
		_drawTertiary = true;
		_border = new Dimension(0, 0);
		_sseq = "";// =
		// "CAGCACGACACUAGCAGUCAGUGUCAGACUGCAIACAGCACGACACUAGCAGUCAGUGUCAGACUGCAIACAGCACGACACUAGCAGUCAGUGUCAGACUGCAIA";
		_sstruct = "";// =
		// "..(((((...(((((...(((((...(((((.....)))))...))))).....(((((...(((((.....)))))...))))).....)))))...)))))..";
		_periodResNum = VARNAConfig.DEFAULT_PERIOD;
		_rotation = 0.0;
		_title = "";
		_titleSize = VARNAConfig.DEFAULT_TITLE_FONT.getSize();
		
		_backboneColor = VARNAConfig.DEFAULT_BACKBONE_COLOR;
		_drawBackbone = true;
		
		_bondColor = VARNAConfig.DEFAULT_BOND_COLOR;
		_bpStyle = VARNAConfig.DEFAULT_BP_STYLE;
		
		_highlightRegion = "";

		_baseOutlineColor = VARNAConfig.BASE_OUTLINE_COLOR_DEFAULT;
		_baseInnerColor = VARNAConfig.BASE_INNER_COLOR_DEFAULT;
		_baseNumColor = VARNAConfig.BASE_NUMBER_COLOR_DEFAULT;
		_baseNameColor = VARNAConfig.BASE_NAME_COLOR_DEFAULT;

		_titleColor = VARNAConfig.DEFAULT_TITLE_COLOR;
		_warning = false;
		_error = true;
		_modifiable = true;
		_zoom = VARNAConfig.DEFAULT_ZOOM;
		_zoomAmount = VARNAConfig.DEFAULT_AMOUNT;

		_comparisonMode = false;
		_firstSequence = "";
		_firstStructure = "";
		_secondSequence = "";
		_secondStructure = "";

		_gapsColor = VARNAConfig.DEFAULT_DASH_BASE_COLOR;
		_useGapsColor = false;
		_nonStandardColor = VARNAConfig.DEFAULT_SPECIAL_BASE_COLOR;
		_useNonStandardColor = false;

		_useInnerBaseColor = false;
		_useBaseNameColor = false;
		_useBaseNumbersColor = false;
		_useBaseOutlineColor = false;

		_bpIncrement = VARNAConfig.DEFAULT_BP_INCREMENT;

		_URL = "";
		_flatExteriorLoop = false;
	}

	public static Color getSafeColor(String col, Color def) {
		Color result;
		try {
			result = Color.decode(col);

		} catch (NumberFormatException e) {
			try {
				result = Color.getColor(col, def);
			} catch (Exception e2) {
				// Not a valid color
				return def;
			}
		}
		return result;
	}

	public static final String LEONTIS_WESTHOF_BP_STYLE = "lw";
	public static final String SIMPLE_BP_STYLE = "simple";
	public static final String RNAVIZ_BP_STYLE = "rnaviz";
	public static final String NONE_BP_STYLE = "none";

	private VARNAConfig.BP_STYLE getSafeBPStyle(String opt,
			VARNAConfig.BP_STYLE def) {
		if (opt.toLowerCase().equals(LEONTIS_WESTHOF_BP_STYLE)) {
			return VARNAConfig.BP_STYLE.BP_STYLE_LW;
		} else if (opt.toLowerCase().equals(SIMPLE_BP_STYLE)) {
			return VARNAConfig.BP_STYLE.BP_STYLE_SIMPLE;
		} else if (opt.toLowerCase().equals(RNAVIZ_BP_STYLE)) {
			return VARNAConfig.BP_STYLE.BP_STYLE_RNAVIZ;
		} else if (opt.toLowerCase().equals(NONE_BP_STYLE)) {
			return VARNAConfig.BP_STYLE.BP_STYLE_NONE;
		} else {
			return def;
		}
	}

	public static String[][] getParameterInfo() {
		String[][] info = {
				// Parameter Name Kind of Value Description
				{
						algoOpt,
						"String",
						"Drawing algorithm, choosen from ["
								+ VARNAConfigLoader.ALGORITHM_NAVIEW + ","
								+ VARNAConfigLoader.ALGORITHM_LINE + ","
								+ VARNAConfigLoader.ALGORITHM_RADIATE + ","
								+ VARNAConfigLoader.ALGORITHM_CIRCULAR + "]" },
				{ annotationsOpt, "string", "A set of textual annotations" },
				{ backboneColorOpt, "Color", "Backbone color (Ex: #334455)" },
				{ drawBackboneOpt, "boolean", "True if the backbone must be drawn, false otherwise" },
				{ backgroundColorOpt, "Color", "Background color (Ex: #334455)" },
				{ baseInnerColorOpt, "Color",
						"Default value for inner base color (Ex: #334455)" },
				{ baseNameColorOpt, "Color",
						"Residues font color (Ex: #334455)" },
				{ baseNumbersColorOpt, "Color",
						"Base numbers font color (Ex: #334455)" },
				{ baseOutlineColorOpt, "Color",
						"Base outline color (Ex: #334455)" },
				{ bondColorOpt, "Color", "Base pair color (Ex: #334455)" },
				{ borderOpt, "String",
						"Border width and height in pixels (Ex: \"20x40\")" },
				{ sequenceOpt, "String", "Raw RNA sequence" },
				{ structureOpt, "String",
						"RNA structure given in dot bracket notation (DBN)" },
				{ periodNumOpt, "int", "Periodicity of base-numbering" },
				{ titleOpt, "String", "RNA drawing title" },
				{ titleColorOpt, "Color", "Title color (Ex: #334455)" },
				{ titleSizeOpt, "int", "Title font size" },

				{ zoomOpt, "int", "Zoom coefficient" },
				{ zoomAmountOpt, "int", "Zoom increment on user interaction" },
				{ warningOpt, "boolean", "Show warnings" },
				{ errorOpt, "boolean", "Show errors" },
				{ gapsBaseColorOpt, "Color","Define and use custom color for gaps bases in comparison mode" },
				{ nonStandardColorOpt, "Color",
						"Define and use custom color for non-standard bases in comparison mode" },
				{ basesStyleOpt, "String", "Base style declaration" },
				{ applyBasesStyleOpt, "String", "Base style application" },
				{
						bpStyleOpt,
						"String",
						"Look and feel for base pairs drawings, choosen from ["
								+ VARNAConfigLoader.LEONTIS_WESTHOF_BP_STYLE
								+ "," + VARNAConfigLoader.NONE_BP_STYLE + ","
								+ VARNAConfigLoader.SIMPLE_BP_STYLE + ","
								+ VARNAConfigLoader.RNAVIZ_BP_STYLE + "]" },
				{ modifiableOpt, "boolean", "Allows/prohibits modifications" },
				{ comparisonModeOpt, "boolean", "Activates comparison mode" },
				{ firstSequenceForComparisonOpt, "String",
						"In comparison mode, sequence of first RNA" },
				{ firstStructureForComparisonOpt, "String",
						"In comparison mode, structure of first RNA" },
				{ secondSequenceForComparisonOpt, "String",
						"In comparison mode, sequence of second RNA" },
				{ secondStructureForComparisonOpt, "String",
						"In comparison mode, structure of second RNA" },
				{
						auxBPsOpt,
						"String",
						"Adds a list of (possibly non-canonical) base-pairs to those already defined by the main secondary structure (Ex: \"(1,10);(2,11);(3,12)\"). Custom BP styles can be specified (Ex: \"(2,11):thickness=4;(3,12):color=#FF0000\")." },
				{
						rotationOpt,
						"String",
						"Rotates RNA after initial drawing (Ex: '20' for a 20 degree counter-clockwise rotation)" },
				{ drawNCOpt, "boolean",
						"Toggles on/off display of non-canonical base-pairs" },
				{ drawTertiaryOpt, "boolean",
						"Toggles on/off display of tertiary interaction, ie pseudoknots" },
				{ flatExteriorLoopOpt, "boolean",
						"Toggles on/off drawing exterior bases on a straight line" },
				{ bpIncrementOpt, "float",
						"Distance between nested base-pairs (i.e. arcs) in linear representation" }, };
		return info;
	}

	private void retrieveParametersValues() throws ExceptionParameterError,
			ExceptionModeleStyleBaseSyntaxError, ExceptionNonEqualLength,
			IOException, ExceptionFileFormatOrSyntax, ExceptionLoadingFailed {
		_numRows = 1;
		_numColumns = 1;
		_basesStyleList = new ArrayList<ModeleStyleBase>();

		try {
			_numRows = Integer.parseInt(_optionProducer.getParameterValue(
					numRowsOpt, "" + _numRows));
		} catch (NumberFormatException e) {
			throw new ExceptionParameterError(e.getMessage(), "'"
					+ _optionProducer.getParameterValue(numRowsOpt, ""
							+ _numRows)
					+ "' is not a integer value for the number of rows !");
		}
		try {
			_numColumns = Integer.parseInt(_optionProducer.getParameterValue(
					numColumnsOpt, "" + _numColumns));
		} catch (NumberFormatException e) {
			throw new ExceptionParameterError(e.getMessage(), "'"
					+ _optionProducer.getParameterValue(numColumnsOpt, ""
							+ _numColumns)
					+ "' is not a integer value for the number of columns !");
		}

		String tmp = null;
		for (int i = 0; i < MAXSTYLE; i++) {
			tmp = _optionProducer.getParameterValue(basesStyleOpt + i, null);
			if (tmp != null) {
				ModeleStyleBase msb = new ModeleStyleBase();
				msb.assignParameters(tmp);
				_basesStyleList.add(msb);
			} else {
				_basesStyleList.add(null);
			}
		}

		// _containerApplet.getLayout().
		int x;
		String n;
		initValues();
		for (int i = 0; i < _numColumns; i++) {
			for (int j = 0; j < _numRows; j++) {
				// initValues();

				x = 1 + j + i * _numRows;
				n = "" + x;
				if ((_numColumns == 1) && (_numRows == 1)) {
					n = "";
				}
				_useGapsColor = false;
				_useNonStandardColor = false;

				_algo = _optionProducer.getParameterValue(algoOpt + n, _algo);
				_customBases = _optionProducer.getParameterValue(customBasesOpt
						+ n, _customBases);
				_auxBPs = _optionProducer.getParameterValue(auxBPsOpt + n,
						_auxBPs);
				_highlightRegion = _optionProducer.getParameterValue(highlightRegionOpt + n, 
						_highlightRegion);

				_chemProbs = _optionProducer.getParameterValue(chemProbOpt + n,
						_chemProbs);
				_customBPs = _optionProducer.getParameterValue(
						customBPsOpt + n, _customBPs);
				_colorMapStyle = _optionProducer.getParameterValue(
						colorMapDefOpt + n, _colorMapStyle);
				_colorMapCaption = _optionProducer.getParameterValue(
						colorMapCaptionOpt + n, _colorMapCaption);
				_backboneColor = getSafeColor(_optionProducer
						.getParameterValue(backboneColorOpt + n, _backboneColor
								.toString()), _backboneColor);
				_backgroundColor = getSafeColor(_optionProducer
						.getParameterValue(backgroundColorOpt + n,
								_backgroundColor.toString()), _backgroundColor);
				_bondColor = getSafeColor(_optionProducer.getParameterValue(
						bondColorOpt + n, _bondColor.toString()), _bondColor);
				_bpStyle = getSafeBPStyle(_optionProducer.getParameterValue(
						bpStyleOpt + n, ""), _bpStyle);
				tmp = _optionProducer.getParameterValue(baseNameColorOpt + n,
						"");
				if (!tmp.equals("")) {
					_useBaseNameColor = true;
					_baseNameColor = getSafeColor(tmp, _baseNameColor);
				}
				tmp = _optionProducer.getParameterValue(
						baseNumbersColorOpt + n, "");
				if (!tmp.equals("")) {
					_useBaseNumbersColor = true;
					_baseNumColor = getSafeColor(tmp, _baseNumColor);
				}
				tmp = _optionProducer.getParameterValue(
						baseOutlineColorOpt + n, "");
				if (!tmp.equals("")) {
					_useBaseOutlineColor = true;
					_baseOutlineColor = getSafeColor(tmp, _baseOutlineColor);
				}
				tmp = _optionProducer.getParameterValue(baseInnerColorOpt + n,
						"");
				if (!tmp.equals("")) {
					_useInnerBaseColor = true;
					_baseInnerColor = getSafeColor(tmp, _baseInnerColor);
				}

				tmp = _optionProducer.getParameterValue(
						nonStandardColorOpt + n, "");
				if (!tmp.equals("")) {
					_nonStandardColor = getSafeColor(tmp, _nonStandardColor);
					_useNonStandardColor = true;
				}
				tmp = _optionProducer.getParameterValue(gapsBaseColorOpt + n,
						_gapsColor.toString());
				if (!tmp.equals("")) {
					_gapsColor = getSafeColor(tmp, _gapsColor);
					_useGapsColor = true;
				}
				try {
					_rotation = Double.parseDouble(_optionProducer
							.getParameterValue(rotationOpt + n, Double
									.toString(_rotation)));
				} catch (NumberFormatException e) {
					throw new ExceptionParameterError(e.getMessage(), "'"
							+ _optionProducer.getParameterValue(periodNumOpt
									+ n, "" + _periodResNum)
							+ "' is not a valid float value for rotation!");
				}

				try {
					_bpIncrement = Double.parseDouble(_optionProducer
							.getParameterValue(bpIncrementOpt + n, Double
									.toString(_bpIncrement)));
				} catch (NumberFormatException e) {
				}

				try {
					_periodResNum = Integer.parseInt(_optionProducer
							.getParameterValue(periodNumOpt + n, ""
									+ _periodResNum));
				} catch (NumberFormatException e) {
					throw new ExceptionParameterError(
							e.getMessage(),
							"'"
									+ _optionProducer.getParameterValue(
											periodNumOpt + n, ""
													+ _periodResNum)
									+ "' is not a valid integer value for the period of residue numbers!");
				}
				_title = _optionProducer
						.getParameterValue(titleOpt + n, _title);
				_titleColor = getSafeColor(_optionProducer.getParameterValue(
						titleColorOpt + n, _titleColor.toString()), _titleColor);
				try {
					_titleSize = Integer.parseInt(_optionProducer
							.getParameterValue(titleSizeOpt + n, ""
									+ _titleSize));
				} catch (NumberFormatException e) {
					throw new ExceptionParameterError(
							e.getMessage(),
							"'"
									+ _optionProducer.getParameterValue(
											titleSizeOpt + n, "" + _titleSize)
									+ "' is not a valid integer value for the number of rows !");
				}

				try {
					_zoom = Double.parseDouble(_optionProducer
							.getParameterValue(zoomOpt + n, "" + _zoom));
				} catch (NumberFormatException e) {
					throw new ExceptionParameterError(e.getMessage(), "'"
							+ _optionProducer.getParameterValue(zoomOpt + n, ""
									+ _zoom)
							+ "' is not a valid integer value for the zoom !");
				}

				try {
					_zoomAmount = Double.parseDouble(_optionProducer
							.getParameterValue(zoomAmountOpt + n, ""
									+ _zoomAmount));
				} catch (NumberFormatException e) {
					throw new ExceptionParameterError(
							e.getMessage(),
							"'"
									+ _optionProducer
											.getParameterValue(zoomAmountOpt
													+ n, "" + _zoomAmount)
									+ "' is not a valid integer value for the zoom amount !");
				}

				_autoHelices = Boolean.parseBoolean(_optionProducer
						.getParameterValue(autoHelicesOpt + n, ""
								+ _autoHelices));

				_drawColorMap = Boolean.parseBoolean(_optionProducer
						.getParameterValue(drawColorMapOpt + n, ""
								+ _drawColorMap));

				_drawBackbone = Boolean.parseBoolean(_optionProducer
						.getParameterValue(drawBackboneOpt + n, ""
								+ _drawBackbone));
				
				_colorMapValues = _optionProducer
				.getParameterValue(colorMapOpt + n, _colorMapValues);

				_autoTerminalLoops = Boolean.parseBoolean(_optionProducer
						.getParameterValue(autoTerminalLoopsOpt + n, ""
								+ _autoTerminalLoops));

				_autoInteriorLoops = Boolean.parseBoolean(_optionProducer
						.getParameterValue(autoInteriorLoopsOpt + n, ""
								+ _autoInteriorLoops));

				_drawNC = Boolean.parseBoolean(_optionProducer
						.getParameterValue(drawNCOpt + n, "" + _drawNC));

				_flatExteriorLoop = Boolean.parseBoolean(_optionProducer
						.getParameterValue(flatExteriorLoopOpt + n, ""
								+ _flatExteriorLoop));

				_drawTertiary = Boolean.parseBoolean(_optionProducer
						.getParameterValue(drawTertiaryOpt + n, ""
								+ _drawTertiary));

				_warning = Boolean.parseBoolean(_optionProducer
						.getParameterValue(warningOpt + n, "false"));

				_error = Boolean.parseBoolean(_optionProducer
						.getParameterValue(errorOpt + n, "true"));

				_border = parseDimension(_optionProducer.getParameterValue(
						borderOpt + n, "0X0"));

				_modifiable = Boolean.parseBoolean(_optionProducer
						.getParameterValue(modifiableOpt + n, "true"));

				_comparisonMode = Boolean.parseBoolean(_optionProducer
						.getParameterValue(comparisonModeOpt + n, "false"));
				// On récupère les deux séquences et structures d'alignement
				_firstSequence = _optionProducer.getParameterValue(
						firstSequenceForComparisonOpt + n, _firstSequence);
				_firstStructure = _optionProducer.getParameterValue(
						firstStructureForComparisonOpt + n, _firstStructure);
				_secondSequence = _optionProducer.getParameterValue(
						secondSequenceForComparisonOpt + n, _secondSequence);
				_secondStructure = _optionProducer.getParameterValue(
						secondStructureForComparisonOpt + n, _secondStructure);
				_annotations = _optionProducer.getParameterValue(annotationsOpt
						+ n, _annotations);
				_URL = _optionProducer.getParameterValue(URLOpt + n, _URL);
				if (!_URL.equals("")) {
					_sstruct = "";
					_sseq = "";
				}
				if (_comparisonMode && _firstSequence != null
						&& _firstStructure != null && _secondSequence != null
						&& _secondStructure != null) {
				} else {
					_sseq = _optionProducer.getParameterValue(sequenceOpt + n,
							_sseq);
					_sstruct = _optionProducer.getParameterValue(structureOpt
							+ n, _sstruct);
					if (!_sseq.equals("") && !_sstruct.equals("")) {
						_URL = "";
					}
					_comparisonMode = false;
				}

				// applique les valeurs des parametres recuperees
				applyValues(n);
			}// fin de boucle sur les lignes
		}// fin de boucle sur les colonnes
	}
	
	private RNA _defaultRNA = new RNA(); 
	
	public void setRNA(RNA r)
	{
		_defaultRNA = r;
	}

	public static final String ALGORITHM_CIRCULAR = "circular";
	public static final String ALGORITHM_NAVIEW = "naview";
	public static final String ALGORITHM_LINE = "line";
	public static final String ALGORITHM_RADIATE = "radiate";
	public static final String ALGORITHM_VARNA_VIEW = "varnaview";
	public static final String ALGORITHM_MOTIF_VIEW = "motifview";

	private void applyValues(String n)
			throws ExceptionParameterError, ExceptionNonEqualLength,
			IOException, ExceptionFileFormatOrSyntax, ExceptionLoadingFailed {
		int algoCode;
		if (_algo.equals(ALGORITHM_CIRCULAR))
			algoCode = RNA.DRAW_MODE_CIRCULAR;
		else if (_algo.equals(ALGORITHM_NAVIEW))
			algoCode = RNA.DRAW_MODE_NAVIEW;
		else if (_algo.equals(ALGORITHM_LINE))
			algoCode = RNA.DRAW_MODE_LINEAR;
		else if (_algo.equals(ALGORITHM_RADIATE))
			algoCode = RNA.DRAW_MODE_RADIATE;
		else if (_algo.equals(ALGORITHM_VARNA_VIEW))
			algoCode = RNA.DRAW_MODE_VARNA_VIEW;
		else if (_algo.equals(ALGORITHM_MOTIF_VIEW))
			algoCode = RNA.DRAW_MODE_MOTIFVIEW;
		else
			algoCode = RNA.DRAW_MODE_RADIATE;

		if (_comparisonMode) {
			_mainSurface = new VARNAPanel(_firstSequence, _firstStructure,
					_secondSequence, _secondStructure, algoCode, "");
		} else {
			if (_URL.equals("")) {
				if (!_sstruct.equals(""))
				{ _mainSurface = new VARNAPanel(_sseq, _sstruct, algoCode); }
				else
				{
					_mainSurface = new VARNAPanel();
					try {
						_defaultRNA.drawRNA(algoCode);
					} catch (ExceptionDrawingAlgorithm e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
					_mainSurface.drawRNA(_defaultRNA);
				}
			} else {
				URL url = null;
				try {
					url = new URL(_URL);
					URLConnection connexion = url.openConnection();
					connexion.setUseCaches(false);
					InputStream r = connexion.getInputStream();
					InputStreamReader inr = new InputStreamReader(r);
					
					RNA rna = new RNA();
					_mainSurface = new VARNAPanel();
					if (_URL.toLowerCase().endsWith(VARNAPanel.VARNA_SESSION_EXTENSION))
					{
						FullBackup f;
							f = VARNAPanel.importSession(r,_URL);
						_mainSurface.setConfig(f.config);
						_mainSurface.showRNA(f.rna);
					}
					else
					{
					rna.loadSecStr(inr, RNA.guessFileTypeFromExtension(_URL));
					rna.drawRNA(algoCode);
					_mainSurface.drawRNA(rna, algoCode);
					}
				} catch (ExceptionFileFormatOrSyntax e) {
					if (url != null)
						e.setPath(url.getPath());
				} catch (ExceptionDrawingAlgorithm e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				if (_mainSurface == null) {
					_mainSurface = new VARNAPanel("X", ".", algoCode);
				}
			}
		}

		if (_useInnerBaseColor) {
			_mainSurface.setBaseInnerColor(_baseInnerColor);
		}
		if (_useBaseOutlineColor) {
			_mainSurface.setBaseOutlineColor(_baseOutlineColor);
		}
		if (_useBaseNameColor) {
			_mainSurface.setBaseNameColor(_baseNameColor);
		}
		if (_useBaseNumbersColor) {
			_mainSurface.setBaseNumbersColor(_baseNumColor);
		}

		_mainSurface.setBackground(_backgroundColor);
		_mainSurface.setNumPeriod(_periodResNum);
		_mainSurface.setBackboneColor(_backboneColor);
		_mainSurface.setDefaultBPColor(_bondColor);
		_mainSurface.setBPHeightIncrement(_bpIncrement);
		_mainSurface.setBPStyle(_bpStyle);
		_mainSurface.setDrawBackbone(_drawBackbone);

		_mainSurface.setTitleFontColor(_titleColor);
		_mainSurface.setTitle(_title);
		_mainSurface.setTitleFontSize(_titleSize);

		_mainSurface.getPopupMenu().get_itemShowWarnings().setState(_warning);
		_mainSurface.setErrorsOn(_error);
		_mainSurface.setFlatExteriorLoop(_flatExteriorLoop);
		_mainSurface.setZoom(_zoom);
		_mainSurface.setZoomIncrement(_zoomAmount);
		_mainSurface.setBorderSize(_border);

		if (_useGapsColor) {
			_mainSurface.setGapsBasesColor(this._gapsColor);
			_mainSurface.setColorGapsBases(true);
		}

		if (_useNonStandardColor) {
			_mainSurface.setNonStandardBasesColor(_nonStandardColor);
			_mainSurface.setColorNonStandardBases(true);
		}

		_mainSurface.setShowNonPlanarBP(_drawTertiary);
		_mainSurface.setShowNonCanonicalBP(_drawNC);
		
		
		applyBasesStyle(n);

		if (!_customBases.equals(""))
			applyBasesCustomStyles(_mainSurface);

		if (!_highlightRegion.equals(""))
			applyHighlightRegion(_mainSurface);

		if (!_auxBPs.equals(""))
			applyAuxBPs(_mainSurface);

		if (!_chemProbs.equals(""))
			applyChemProbs(_mainSurface);

		if (!_customBPs.equals(""))
			applyBPsCustomStyles(_mainSurface);

		_mainSurface.drawRNA();

		if (!_annotations.equals(""))
			applyAnnotations(_mainSurface);
		if (_autoHelices)
			_mainSurface.getVARNAUI().UIAutoAnnotateHelices();
		if (_autoTerminalLoops)
			_mainSurface.getVARNAUI().UIAutoAnnotateTerminalLoops();
		if (_autoInteriorLoops)
			_mainSurface.getVARNAUI().UIAutoAnnotateInteriorLoops();
		_mainSurface.globalRotation(_rotation);

		_mainSurface.setModifiable(_modifiable);
		
		_mainSurface.setColorMapCaption(_colorMapCaption);
		applyColorMapStyle(_mainSurface);
		applyColorMapValues(_mainSurface);
		//if (!_drawColorMap)
		//  _mainSurface.drawColorMap(_drawColorMap);

		// ajoute le VARNAPanel au conteneur
		_VARNAPanelList.add(_mainSurface);
	}
	
	private void applyBasesStyle(String n) throws ExceptionParameterError
	{
		String tmp;
		for (int numStyle = 0; numStyle < _basesStyleList.size(); numStyle++) 
		{
			if (_basesStyleList.get(numStyle) != null) {
				tmp = _optionProducer.getParameterValue(applyBasesStyleOpt + numStyle + "on" + n, null);
				ArrayList<Integer> indicesList = new ArrayList<Integer>();
				if (tmp != null) {
					String[] basesList = tmp.split(",");
					for (int k = 0; k < basesList.length; k++) 
					{
						String cand = basesList[k].trim();
						try {
							String[] args = cand.split("-");
							if (args.length==1)
							{
								int baseNum = Integer.parseInt(cand);
								int index = _mainSurface.getRNA().getIndexFromBaseNumber(baseNum);
								if (index!=-1)
								{ indicesList.add(index); }
							}
							else if (args.length==2)
							{
								int baseNumFrom = Integer.parseInt(args[0].trim());
								int indexFrom = _mainSurface.getRNA().getIndexFromBaseNumber(baseNumFrom);
								int baseNumTo = Integer.parseInt(args[1].trim());
								int indexTo = _mainSurface.getRNA().getIndexFromBaseNumber(baseNumTo);
								if ((indexFrom!=-1)&&(indexTo!=-1))
								{ 
									for (int l=indexFrom;l<=indexTo;l++)
										indicesList.add(l); 
								}								
							}
						} catch (NumberFormatException e) {
							throw new ExceptionParameterError(e
									.getMessage(), "Bad Base Index: "
									+ basesList[k]);
						}				
					}
					for(int k=0;k<indicesList.size();k++)
					{
						int index = indicesList.get(k);
					if ((index >= 0) && (index < _mainSurface.getRNA().get_listeBases().size())) 
					{
						_mainSurface.getRNA().get_listeBases().get(
								index).setStyleBase(
								_basesStyleList.get(numStyle));
					}
					}
				}
			}
		}// fin de boucle sur les styles

	}
	
	private void applyColorMapStyle(VARNAPanel vp)
	{
		if (_colorMapStyle.length()==0)
			return;
		vp.setColorMap(ModeleColorMap.parseColorMap(_colorMapStyle));
	}
	
	private void applyColorMapValues(VARNAPanel vp)
	{
		if (!_colorMapValues.equals(""))
		{	  
		String[] values = _colorMapValues.split("[;,]");
		ArrayList<Double> vals = new ArrayList<Double>();
		for (int i = 0; i < values.length; i++) {
			try{
			vals.add(Double.parseDouble(values[i]));
			}
			catch(Exception e)
			{
				
			}
		}
		Double[] result = new Double[vals.size()];
		vals.toArray(result);
		vp.setColorMapValues(result);
		_drawColorMap = true;
		}
	}
	

	private void applyBasesCustomStyles(VARNAPanel vp) {
		String[] baseStyles = _customBases.split(";");
		for (int i = 0; i < baseStyles.length; i++) {
			String thisStyle = baseStyles[i];
			String[] data = thisStyle.split(":");
			try {
				if (data.length == 2) {
					int baseNum = Integer.parseInt(data[0]);
					int index = _mainSurface.getRNA().getIndexFromBaseNumber(baseNum);
					if (index!=-1)
					{
						String style = data[1];
						ModeleStyleBase msb = vp.getRNA().get_listeBases().get(index).getStyleBase().clone();
						msb.assignParameters(style);
						vp.getRNA().get_listeBases().get(index).setStyleBase(msb);
					}
				}
			} catch (Exception e) {
				System.err.println("ApplyBasesCustomStyle: " + e.toString());
			}
		}
	}
	
	private void applyHighlightRegion(VARNAPanel vp) {
		String[] regions = _highlightRegion.split(";");
		for (int i = 0; i < regions.length; i++) 
		{
			String region = regions[i];
			System.out.println(region);
			try{
				HighlightRegionAnnotation nt = HighlightRegionAnnotation.parseHighlightRegionAnnotation(region,vp);
				if (nt!=null)
					vp.addHighlightRegion(nt);
			} catch (Exception e) {
				System.err.println("applyHighlightRegion: " + e.toString());
			}
		}
	}
	
	

	private Dimension parseDimension(String s) {
		Dimension d = new Dimension(0, 0);
		try {
			s = s.toLowerCase();
			int i = s.indexOf('x');
			String w = s.substring(0, i);
			String h = s.substring(i + 1);
			d.width = Integer.parseInt(w);
			d.height = Integer.parseInt(h);
		} catch (NumberFormatException e) {
		}
		return d;
	}

	private void applyBPsCustomStyles(VARNAPanel vp) {
		String[] baseStyles = _customBPs.split(";");
		for (int i = 0; i < baseStyles.length; i++) {
			String thisStyle = baseStyles[i];
			String[] data = thisStyle.split(":");
			try {
				if (data.length == 2) {
					String indices = data[0];
					String style = data[1];
					String[] data2 = indices.split(",");
					if (data2.length == 2) {
						String s1 = data2[0];
						String s2 = data2[1];
						if (s1.startsWith("(") && s2.endsWith(")")) {
							int a = Integer.parseInt(s1.substring(1));
							int b = Integer.parseInt(s2.substring(0, s2
									.length() - 1));
							ModeleStyleBP msbp = vp.getRNA().getBPStyle(a, b);
							if (msbp != null) {
								msbp.assignParameters(style);
							}
						}
					}
				}
			} catch (Exception e) {
				System.err.println("ApplyBPsCustomStyle: " + e.toString());
			}
		}
	}

	private void applyChemProbs(VARNAPanel vp) {
		String[] chemProbs = _chemProbs.split(";");
		for (int i = 0; i < chemProbs.length; i++) {
			String thisAnn = chemProbs[i];
			String[] data = thisAnn.split(":");
			try {
				if (data.length == 2) {
					String indices = data[0];
					String style = data[1];
					String[] data2 = indices.split("-");
					if (data2.length == 2) {
						int a = Integer.parseInt(data2[0]);
						int b = Integer.parseInt(data2[1]);
						int c = vp.getRNA().getIndexFromBaseNumber(a);
						int d = vp.getRNA().getIndexFromBaseNumber(b);
						ArrayList<ModeleBase> mbl = vp.getRNA()
								.get_listeBases();
						ChemProbAnnotation cpa = new ChemProbAnnotation(mbl
								.get(c), mbl.get(d), style);
						vp.getRNA().addChemProbAnnotation(cpa);
					}
				}
			} catch (Exception e) {
				System.err.println("ChempProbs: " + e.toString());
			}
		}
	}

	private void applyAuxBPs(VARNAPanel vp) {
		String[] baseStyles = _auxBPs.split(";");
		for (int i = 0; i < baseStyles.length; i++) {
			String thisStyle = baseStyles[i];
			String[] data = thisStyle.split(":");
			try {
				if (data.length >= 1) {
					String indices = data[0];
					String[] data2 = indices.split(",");
					if (data2.length == 2) {
						String s1 = data2[0];
						String s2 = data2[1];
						if (s1.startsWith("(") && s2.endsWith(")")) {
							int a = Integer.parseInt(s1.substring(1));
							int b = Integer.parseInt(s2.substring(0, s2
									.length() - 1));
							int c = vp.getRNA().getIndexFromBaseNumber(a);
							int d = vp.getRNA().getIndexFromBaseNumber(b);
							ModeleStyleBP msbp = new ModeleStyleBP(vp.getRNA()
									.get_listeBases().get(c), vp.getRNA()
									.get_listeBases().get(d));
							if (data.length >= 2) {
								String style = data[1];
								msbp.assignParameters(style);
							}
							vp.getRNA().addBPToStructure(a, b, msbp);
						}
					}
				}
			} catch (Exception e) {
				System.err.println("AuxApplyBPs: " + e.toString());
			}
		}
	}

	/**
	 * Format:
	 * string:[type=[H|B|L|P]|x=double|y=double|anchor=int|size=int|color=Color];
	 * 
	 * @param vp
	 */
	private void applyAnnotations(VARNAPanel vp) {
		String[] annotations = _annotations.split(";");
		for (int i = 0; i < annotations.length; i++) {
			String thisAnn = annotations[i];
			String[] data = thisAnn.split(":");

			String text = "";
			int anchor = -1;
			int x = -1;
			int y = -1;
			int type = TextAnnotation.LOOP;
			Font font = TextAnnotation.DEFAULTFONT;
			Color color = TextAnnotation.DEFAULTCOLOR;
			TextAnnotation ann = null;
			try {
				if (data.length == 2) {
					text = data[0];
					String[] data2 = data[1].split(",");
					for (int j = 0; j < data2.length; j++) {
						String opt = data2[j];
						String[] data3 = opt.split("=");
						if (data3.length == 2) {
							String name = data3[0].toLowerCase();
							String value = data3[1];
							if (name.equals("type")) {
								if (value.toUpperCase().equals("H")) {
									type = TextAnnotation.HELIX;
								} else if (value.toUpperCase().equals("L")) {
									type = TextAnnotation.LOOP;
								} else if (value.toUpperCase().equals("P")) {
									type = TextAnnotation.POSITION;
								} else if (value.toUpperCase().equals("B")) {
									type = TextAnnotation.BASE;
								}
							} else if (name.equals("x")) {
								x = Integer.parseInt(value);
							} else if (name.equals("y")) {
								y = Integer.parseInt(value);
							} else if (name.equals("anchor")) {
								anchor = Integer.parseInt(value);
							} else if (name.equals("size")) {
								font = font.deriveFont((float) Integer
										.parseInt(value));
							} else if (name.equals("color")) {
								color = getSafeColor(value, color);
							}
						}
					}
					switch (type) {
					case (TextAnnotation.POSITION):
						if ((x != -1) && (y != -1)) {
							Point2D.Double p = vp
									.panelToLogicPoint(new Point2D.Double(x, y));
							ann = new TextAnnotation(text, p.x, p.y);
						}
						break;
					case (TextAnnotation.BASE):
						if (anchor != -1) {
							int index = vp.getRNA().getIndexFromBaseNumber(
									anchor);
							ModeleBase mb = vp.getRNA().get_listeBases().get(
									index);
							ann = new TextAnnotation(text, mb);
						}
						break;
					case (TextAnnotation.HELIX):
						if (anchor != -1) {
							ArrayList<ModeleBase> mbl = new ArrayList<ModeleBase>();
							int index = vp.getRNA().getIndexFromBaseNumber(
									anchor);
							ArrayList<Integer> il = vp.getRNA()
									.findHelix(index);
							for (int k : il) {
								mbl.add(vp.getRNA().get_listeBases().get(k));
							}
							ann = new TextAnnotation(text, mbl, type);
						}
						break;
					case (TextAnnotation.LOOP):
						if (anchor != -1) {
							ArrayList<ModeleBase> mbl = new ArrayList<ModeleBase>();
							int index = vp.getRNA().getIndexFromBaseNumber(
									anchor);
							ArrayList<Integer> il = vp.getRNA().findLoop(index);
							for (int k : il) {
								mbl.add(vp.getRNA().get_listeBases().get(k));
							}
							ann = new TextAnnotation(text, mbl, type);
						}
						break;
					}
					if (ann != null) {
						ann.setColor(color);
						ann.setFont(font);
						vp.addAnnotation(ann);
					}
				}
			} catch (Exception e) {
				System.err.println("Apply Annotations: " + e.toString());
			}
		}
	}
}
