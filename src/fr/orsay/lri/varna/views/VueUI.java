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
package fr.orsay.lri.varna.views;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics2D;
import java.awt.Image;
import java.awt.Point;
import java.awt.Toolkit;
import java.awt.geom.AffineTransform;
import java.awt.geom.Point2D;
import java.awt.image.BufferedImage;
import java.awt.image.FilteredImageSource;
import java.awt.image.ImageFilter;
import java.awt.image.ImageProducer;
import java.awt.image.RGBImageFilter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Vector;

import javax.imageio.IIOImage;
import javax.imageio.ImageIO;
import javax.imageio.ImageWriteParam;
import javax.imageio.ImageWriter;
import javax.imageio.stream.FileImageOutputStream;
import javax.swing.JColorChooser;
import javax.swing.JFileChooser;
import javax.swing.JOptionPane;
import javax.swing.filechooser.FileFilter;

import fr.orsay.lri.varna.VARNAPanel;
import fr.orsay.lri.varna.applications.FileNameExtensionFilter;
import fr.orsay.lri.varna.applications.VARNAPrinter;
import fr.orsay.lri.varna.exceptions.ExceptionExportFailed;
import fr.orsay.lri.varna.exceptions.ExceptionFileFormatOrSyntax;
import fr.orsay.lri.varna.exceptions.ExceptionJPEGEncoding;
import fr.orsay.lri.varna.exceptions.ExceptionLoadingFailed;
import fr.orsay.lri.varna.exceptions.ExceptionNonEqualLength;
import fr.orsay.lri.varna.exceptions.ExceptionPermissionDenied;
import fr.orsay.lri.varna.exceptions.ExceptionUnmatchedClosingParentheses;
import fr.orsay.lri.varna.exceptions.ExceptionWritingForbidden;
import fr.orsay.lri.varna.models.VARNAConfig;
import fr.orsay.lri.varna.models.annotations.ChemProbAnnotation;
import fr.orsay.lri.varna.models.annotations.HighlightRegionAnnotation;
import fr.orsay.lri.varna.models.annotations.TextAnnotation;
import fr.orsay.lri.varna.models.rna.ModeleBase;
import fr.orsay.lri.varna.models.rna.ModeleBaseNucleotide;
import fr.orsay.lri.varna.models.rna.ModeleBasesComparison;
import fr.orsay.lri.varna.models.rna.ModeleStyleBP;
import fr.orsay.lri.varna.models.rna.RNA;

public class VueUI {
	private VARNAPanel _vp;
	private File _fileChooserDirectory = null;

	public VueUI(VARNAPanel vp) {
		_vp = vp;
	}

	public void UIRadiate() {
		if (_vp.isModifiable()) {
			_vp.reset();
			_vp.drawRNA(_vp.getRNA(), RNA.DRAW_MODE_RADIATE);
			_vp.repaint();
		}
	}

	public void UIToggleColorMap() {
		if (_vp.isModifiable()) {
			_vp.setColorMapVisible(!_vp.getColorMapVisible());
			_vp.repaint();
		}
	}
	
	public void UIToggleFlatExteriorLoop() {
		if (_vp.isModifiable()
				&& _vp.getRNA().get_drawMode() == RNA.DRAW_MODE_RADIATE) {
			_vp.setFlatExteriorLoop(!_vp.getFlatExteriorLoop());
			UIRadiate();
		}
	}
	
	public void UIMOTIFView() {
		if (_vp.isModifiable()) {
			_vp.reset();
			_vp.drawRNA(_vp.getRNA(), RNA.DRAW_MODE_MOTIFVIEW);
			_vp.repaint();
		}
	}

	public void UILine() {
		if (_vp.isModifiable()) {
			_vp.reset();
			_vp.drawRNA(_vp.getRNA(), RNA.DRAW_MODE_LINEAR);
			_vp.repaint();
		}
	}

	public void UICircular() {
		if (_vp.isModifiable()) {
			_vp.reset();
			_vp.drawRNA(_vp.getRNA(), RNA.DRAW_MODE_CIRCULAR);
			_vp.repaint();
		}
	}

	public void UINAView() {
		if (_vp.isModifiable()) {
			_vp.reset();
			_vp.drawRNA(_vp.getRNA(), RNA.DRAW_MODE_NAVIEW);
			_vp.repaint();
		}
	}

	public void UIVARNAView() {
		if (_vp.isModifiable()) {
			_vp.reset();
			_vp.drawRNA(_vp.getRNA(), RNA.DRAW_MODE_VARNA_VIEW);
			_vp.repaint();
		}
	}

	public void UIReset() {
		if (_vp.isModifiable()) {
			_vp.reset();
			_vp.drawRNA(_vp.getRNA(), _vp.getRNA().get_drawMode());
			_vp.repaint();
		}
	}
	
	private void savePath(JFileChooser jfc)
	{
		_fileChooserDirectory = jfc.getCurrentDirectory();
	}

	private void loadPath(JFileChooser jfc)
	{
		if (_fileChooserDirectory != null)
		{
			jfc.setCurrentDirectory(_fileChooserDirectory);
		}
	}
	
	
	public void UIFile() throws ExceptionNonEqualLength {
		if (_vp.isModifiable()) {
			JFileChooser fc = new JFileChooser();
			fc.setFileSelectionMode(JFileChooser.OPEN_DIALOG);
			fc.setDialogTitle("Open...");
			loadPath(fc);
			if (fc.showOpenDialog(_vp) == JFileChooser.APPROVE_OPTION) {
				try {
					savePath(fc);
					String path = fc.getSelectedFile().getAbsolutePath();
					if (!path.toLowerCase().endsWith(".varna"))
					{
						RNA r = new RNA();
						r.loadSecStr(path);
						_vp.drawRNAInterpolated(r);
						_vp.repaint();
					}
					else
					{
						_vp.loadSession(path);
					}
				} catch (ExceptionExportFailed e1) {
					_vp.errorDialog(e1);
				} catch (ExceptionPermissionDenied e1) {
					_vp.errorDialog(e1);
				} catch (ExceptionLoadingFailed e1) {
					_vp.errorDialog(e1);
				} catch (ExceptionFileFormatOrSyntax e1) {
					_vp.errorDialog(e1);
				} catch (ExceptionUnmatchedClosingParentheses e1) {
					_vp.errorDialog(e1);
				} catch (FileNotFoundException e) {
					_vp.errorDialog(e);
				}
			}
		}
	}

	public void UISetColorMapStyle() {
			VueColorMapStyle cms = new VueColorMapStyle(_vp);
			if (JOptionPane.showConfirmDialog(_vp, cms,
					"Choose color map style", JOptionPane.OK_CANCEL_OPTION) == JOptionPane.OK_OPTION) {
				_vp.setColorMap(cms.getColorMap());
			}
			else
			{
				cms.cancelChanges();
			}
	}

	public void UILoadColorMapValues() {
		VueLoadColorMapValues cms = new VueLoadColorMapValues(_vp);
		if (JOptionPane.showConfirmDialog(_vp, cms,
				"Load base values", JOptionPane.OK_CANCEL_OPTION) == JOptionPane.OK_OPTION) {
			try {
				_vp.setColorMapVisible(true);
				_vp.readValues(cms.getReader());
			} catch (IOException e) {
				_vp.errorDialog(e);
			}

		}
	}
	
	
	public void UISetColorMapValues() {
		VueBaseValues cms = new VueBaseValues(_vp);
		if (JOptionPane.showConfirmDialog(_vp, cms,
				"Choose base values", JOptionPane.OK_CANCEL_OPTION) == JOptionPane.OK_OPTION) {	
		}
		else
		{
			cms.cancelChanges();
		}
	}
	
	public void UIManualInput() throws ParseException, ExceptionNonEqualLength {
		if (_vp.isModifiable()) {
			VueManualInput manualInput = new VueManualInput(_vp);
			if (JOptionPane.showConfirmDialog(_vp, manualInput.getPanel(),
					"Input sequence/structure", JOptionPane.OK_CANCEL_OPTION) == JOptionPane.OK_OPTION) {
				if (_vp.getRNA().getSize() == 0) {

				}
				_vp.drawRNAInterpolated(manualInput.getTseq().getText(),
						manualInput.getTstr().getText(), _vp.getRNA()
								.get_drawMode());
				_vp.repaint();
			}
		}
	}

	public void UISetTitle() {
		if (_vp.isModifiable()) {
			String res = JOptionPane.showInputDialog(_vp, "Input title", _vp
					.getTitle());
			if (res != null) {
				_vp.setTitle(res);
				_vp.repaint();
			}
		}
	}

	public void UISetColorMapCaption() {
		if (_vp.isModifiable()) {
			String res = JOptionPane.showInputDialog(_vp, "Input new color map caption", _vp.getColorMapCaption());
			if (res != null) {
				_vp.setColorMapCaption(res);
				_vp.repaint();
			}
		}
	}

	public void UISetBaseCharacter() {
		if (_vp.isModifiable()) {
			if (_vp.isComparisonMode()) {
				String res = JOptionPane.showInputDialog(_vp, "Input base",
						((ModeleBasesComparison) _vp.getRNA().get_listeBases()
								.get(_vp.getNearestBase())).getBases());
				if (res != null) {
					((ModeleBasesComparison) _vp.getRNA().get_listeBases().get(
							_vp.getNearestBase())).set_base1(res.charAt(0));
					((ModeleBasesComparison) _vp.getRNA().get_listeBases().get(
							_vp.getNearestBase())).set_base2(res.charAt(1));
					_vp.repaint();
				}

			} else {
				String res = JOptionPane.showInputDialog(_vp, "Input base",
						((ModeleBaseNucleotide) _vp.getRNA().get_listeBases()
								.get(_vp.getNearestBase())).get_c());
				if (res != null) {
					((ModeleBaseNucleotide) _vp.getRNA().get_listeBases().get(
							_vp.getNearestBase())).set_c(res.charAt(0));
					_vp.repaint();
				}
			}
		}
	}

	FileNameExtensionFilter _varnaFilter = new FileNameExtensionFilter(
			"VARNA Session File", "varna", "VARNA");
	FileNameExtensionFilter _bpseqFilter = new FileNameExtensionFilter(
			"BPSeq (CRW) File", "bpseq", "BPSEQ");
	FileNameExtensionFilter _ctFilter = new FileNameExtensionFilter(
			"Connect (MFold) File", "ct", "CT");
	FileNameExtensionFilter _dbnFilter = new FileNameExtensionFilter(
			"Dot-bracket notation (Vienna) File", "dbn", "DBN", "faa", "FAA");

	FileNameExtensionFilter _jpgFilter = new FileNameExtensionFilter(
			"JPEG Picture", "jpeg", "jpg", "JPG", "JPEG");
	FileNameExtensionFilter _pngFilter = new FileNameExtensionFilter(
			"PNG Picture", "png", "PNG");
	FileNameExtensionFilter _epsFilter = new FileNameExtensionFilter(
			"EPS File", "eps", "EPS");
	FileNameExtensionFilter _svgFilter = new FileNameExtensionFilter(
			"SVG Picture", "svg", "SVG");
	FileNameExtensionFilter _xfigFilter = new FileNameExtensionFilter(
			"XFig Diagram", "fig", "xfig", "FIG", "XFIG");

	public void UIExport() throws ExceptionExportFailed,
			ExceptionPermissionDenied, ExceptionWritingForbidden,
			ExceptionJPEGEncoding {
		ArrayList<FileNameExtensionFilter> v = new ArrayList<FileNameExtensionFilter>();
		v.add(_epsFilter);
		v.add(_svgFilter);
		v.add(_xfigFilter);
		v.add(_jpgFilter);
		v.add(_pngFilter);
		String dest = UIChooseOutputFile(v);
		if (dest != null) {
			String extLower = dest.substring(dest.lastIndexOf('.'))
					.toLowerCase();
			// System.out.println(extLower);
			if (extLower.equals(".eps")) {
				_vp.getRNA().saveRNAEPS(dest, _vp.getConfig());
			} else if (extLower.equals(".svg")) {
				_vp.getRNA().saveRNASVG(dest, _vp.getConfig());
			} else if (extLower.equals(".fig") || extLower.equals("xfig")) {
				_vp.getRNA().saveRNAXFIG(dest, _vp.getConfig());
			} else if (extLower.equals(".png")) {
				saveToPNG(dest);
			} else if (extLower.equals("jpg") || extLower.equals("jpeg")) {
				saveToJPEG(dest);
			}
		}
	}

	public void UIExportJPEG() throws ExceptionJPEGEncoding,
			ExceptionExportFailed {
		String dest = UIChooseOutputFile(_jpgFilter);
		if (dest != null) {
			saveToJPEG(dest);
		}
	}

	public void UIPrint() {
		VARNAPrinter.printComponent(_vp);
	}
	
	public void UIExportPNG() throws ExceptionExportFailed {
		String dest = UIChooseOutputFile(_pngFilter);
		if (dest != null) {
			saveToPNG(dest);
		}
	}

	public void UIExportXFIG() throws ExceptionExportFailed,
			ExceptionWritingForbidden {
		String dest = UIChooseOutputFile(_xfigFilter);
		if (dest != null) {
			_vp.getRNA().saveRNAXFIG(dest, _vp.getConfig());
		}
	}

	public void UIExportEPS() throws ExceptionExportFailed,
			ExceptionWritingForbidden {
		String dest = UIChooseOutputFile(_epsFilter);
		if (dest != null) {
			_vp.getRNA().saveRNAEPS(dest, _vp.getConfig());
		}
	}

	public void UIExportSVG() throws ExceptionExportFailed,
			ExceptionWritingForbidden {
		String dest = UIChooseOutputFile(_svgFilter);
		if (dest != null) {
			_vp.getRNA().saveRNASVG(dest, _vp.getConfig());
		}
	}

	public void UISaveAsDBN() throws ExceptionExportFailed,
			ExceptionPermissionDenied {
		String name = _vp.getVARNAUI().UIChooseOutputFile(_dbnFilter);
		if (name != null)
			_vp.getRNA().saveAsDBN(name, _vp.getTitle());
	}

	public void UISaveAsCT() throws ExceptionExportFailed,
			ExceptionPermissionDenied {
		String name = _vp.getVARNAUI().UIChooseOutputFile(_ctFilter);
		if (name != null)
			_vp.getRNA().saveAsCT(name, _vp.getTitle());
	}

	public void UISaveAsBPSEQ() throws ExceptionExportFailed,
			ExceptionPermissionDenied {
		String name = _vp.getVARNAUI().UIChooseOutputFile(_bpseqFilter);
		if (name != null)
			_vp.getRNA().saveAsBPSEQ(name, _vp.getTitle());
	}

	public void UISaveAs() throws ExceptionExportFailed,
			ExceptionPermissionDenied {
		ArrayList<FileNameExtensionFilter> v = new ArrayList<FileNameExtensionFilter>();
		v.add(_bpseqFilter);
		v.add(_dbnFilter);
		v.add(_ctFilter);
		v.add(_varnaFilter);
		String dest = UIChooseOutputFile(v);
		if (dest != null) {
			String extLower = dest.substring(dest.lastIndexOf('.'))
					.toLowerCase();
			if (extLower.endsWith("bpseq")) {
				_vp.getRNA().saveAsBPSEQ(dest, _vp.getTitle());
			} else if (extLower.endsWith("ct")) {
				_vp.getRNA().saveAsCT(dest, _vp.getTitle());
			} else if (extLower.endsWith("dbn") || extLower.endsWith("faa")) {
				_vp.getRNA().saveAsDBN(dest, _vp.getTitle());
			} else if (extLower.endsWith("varna")) {
				_vp.saveSession(dest);
			}
		}
	}

	public String UIChooseOutputFile(FileNameExtensionFilter filtre) {
		ArrayList<FileNameExtensionFilter> v = new ArrayList<FileNameExtensionFilter>();
		v.add(filtre);
		return UIChooseOutputFile(v);
	}

	/**
	 * Opens a save dialog with right extensions and return the absolute path
	 * 
	 * @param filtre
	 *            Allowed extensions
	 * @return <code>null</code> if the user doesn't approve the save dialog,<br>
	 *         <code>absolutePath</code> if the user approve the save dialog
	 */
	public String UIChooseOutputFile(ArrayList<FileNameExtensionFilter> filtre) {
		JFileChooser fc = new JFileChooser();
		loadPath(fc);
		String absolutePath = null;
		// applique le filtre
		for (int i = 0; i < filtre.size(); i++) {
			fc.addChoosableFileFilter(filtre.get(i));
		}
		// en mode open dialog pour voir les autres fichiers avec la meme
		// extension
		fc.setFileSelectionMode(JFileChooser.OPEN_DIALOG);
		fc.setDialogTitle("Save...");
		// Si l'utilisateur a valider
		if (fc.showSaveDialog(_vp) == JFileChooser.APPROVE_OPTION) {
			savePath(fc);
			absolutePath = fc.getSelectedFile().getAbsolutePath();
			String extension = _vp.getPopupMenu().get_controleurMenu()
					.getExtension(fc.getSelectedFile());
			FileFilter f = fc.getFileFilter();
			if (f instanceof FileNameExtensionFilter) {
				ArrayList<String> listeExtension = new ArrayList<String>();
				listeExtension.addAll(Arrays
						.asList(((FileNameExtensionFilter) f).getExtensions()));
				// si l'extension du fichier ne fait pas partie de la liste
				// d'extensions acceptées
				if (!listeExtension.contains(extension)) {
					absolutePath += "." + listeExtension.get(0);
				}
			}
		}
		return absolutePath;
	}

	public void UISetBorder() {
		VueBorder border = new VueBorder(_vp);
		Dimension oldBorder = _vp.getBorderSize();
		_vp.drawBBox(true);
		_vp.drawBorder(true);
		_vp.repaint();
		if (JOptionPane.showConfirmDialog(_vp, border.getPanel(),
				"Set new border size", JOptionPane.OK_CANCEL_OPTION) != JOptionPane.OK_OPTION) {
			_vp.setBorderSize(oldBorder);
		}
		_vp.drawBorder(false);
		_vp.drawBBox(false);
		_vp.repaint();
	}

	public void UISetBackground() {
		Color c = JColorChooser.showDialog(_vp, "Choose new background color",
				_vp.getBackground());
		if (c != null) {
			_vp.setBackground(c);
			_vp.repaint();
		}
	}

	public void UIZoomIn() {
		double _actualZoom = _vp.getZoom();
		double _actualAmount = _vp.getZoomIncrement();
		Point _actualTranslation = _vp.getTranslation();
		double newZoom = Math.min(VARNAConfig.MAX_ZOOM, _actualZoom
				* _actualAmount);
		double ratio = newZoom / _actualZoom;
		Point newTrans = new Point((int) (_actualTranslation.x * ratio),
				(int) (_actualTranslation.y * ratio));
		_vp.setZoom(newZoom);
		_vp.setTranslation(newTrans);
		// verification que la translation ne pose pas de problemes
		_vp.checkTranslation();
		_vp.repaint();
	}

	public void UIZoomOut() {
		double _actualZoom = _vp.getZoom();
		double _actualAmount = _vp.getZoomIncrement();
		Point _actualTranslation = _vp.getTranslation();
		double newZoom = Math.max(_actualZoom / _actualAmount,
				VARNAConfig.MIN_ZOOM);
		double ratio = newZoom / _actualZoom;
		Point newTrans = new Point((int) (_actualTranslation.x * ratio),
				(int) (_actualTranslation.y * ratio));
		_vp.setTranslation(newTrans);
		_vp.setZoom(newZoom);
		// verification que la translation ne pose pas de problemes
		_vp.checkTranslation();
		_vp.repaint();
	}

	public void UICustomZoom() {
		VueZoom zoom = new VueZoom(_vp);
		double oldZoom = _vp.getZoom();
		double oldZoomAmount = _vp.getZoomIncrement();
		_vp.drawBBox(true);
		_vp.repaint();
		if (JOptionPane.showConfirmDialog(_vp, zoom.getPanel(), "Set zoom",
				JOptionPane.OK_CANCEL_OPTION) != JOptionPane.OK_OPTION) {
			_vp.setZoom(oldZoom);
			_vp.setZoomIncrement(oldZoomAmount);
		}
		_vp.drawBBox(false);
		_vp.repaint();
	}

	public void UIGlobalRotation() {
		if (_vp.getRNA().get_listeBases().size() > 0) {
			_vp.drawBBox(true);
			_vp.repaint();
			VueGlobalRotation rotation = new VueGlobalRotation(_vp);
			if (JOptionPane.showConfirmDialog(_vp, rotation.getPanel(),
					"Rotate the whole RNA", JOptionPane.OK_CANCEL_OPTION) != JOptionPane.OK_OPTION) {
				_vp.globalRotation(-rotation.getAngle());
			}
			_vp.drawBBox(false);
			_vp.repaint();
		}
	}

	public void UISetBPStyle() {
		if (_vp.getRNA().get_listeBases().size() > 0) {
			VueStyleBP bpstyle = new VueStyleBP(_vp);
			VARNAConfig.BP_STYLE bck = _vp.getBPStyle();
			if (JOptionPane.showConfirmDialog(_vp, bpstyle.getPanel(),
					"Set main base pair style", JOptionPane.OK_CANCEL_OPTION) != JOptionPane.OK_OPTION) {
				_vp.setBPStyle(bck);
				_vp.repaint();
			}
		}
	}

	public void UISetTitleColor() {
		if (_vp.isModifiable()) {

			Color c = JColorChooser.showDialog(_vp, "Choose new title color",
					_vp.getTitleColor());
			if (c != null) {
				_vp.setTitleColor(c);
				_vp.repaint();
			}
		}
	}

	public void UISetBackboneColor() {
		if (_vp.isModifiable()) {
			Color c = JColorChooser.showDialog(_vp,
					"Choose new backbone color", _vp.getBackboneColor());
			if (c != null) {
				_vp.setBackboneColor(c);
				_vp.repaint();
			}
		}
	}

	public void UISetTitleFont() {
		if (_vp.isModifiable()) {
			VueFont font = new VueFont(_vp);
			if (JOptionPane.showConfirmDialog(_vp, font.getPanel(),
					"New Title font", JOptionPane.OK_CANCEL_OPTION) == JOptionPane.OK_OPTION) {
				_vp.setTitleFont(font.getFont());
				_vp.repaint();
			}
		}
	}

	public void UISetSpaceBetweenBases() {
		if (_vp.isModifiable()) {

			VueSpaceBetweenBases vsbb = new VueSpaceBetweenBases(_vp);
			Double oldSpace = _vp.getRNA().get_spaceBetweenBases();
			if (JOptionPane.showConfirmDialog(_vp, vsbb.getPanel(),
					"Set the space between each base",
					JOptionPane.OK_CANCEL_OPTION) != JOptionPane.OK_OPTION) {
				_vp.getRNA().set_spaceBetweenBases(oldSpace);
				_vp.drawRNA(_vp.getRNA());
				_vp.repaint();
			}
		}
	}

	public void UISetBPHeightIncrement() {
		if (_vp.isModifiable()) {

			VueBPHeightIncrement vsbb = new VueBPHeightIncrement(_vp);
			Double oldSpace = _vp.getBPHeightIncrement();
			if (JOptionPane.showConfirmDialog(_vp, vsbb.getPanel(),
					"Set the vertical increment in linear mode",
					JOptionPane.OK_CANCEL_OPTION) != JOptionPane.OK_OPTION) {
				_vp.setBPHeightIncrement(oldSpace);
				_vp.drawRNA(_vp.getRNA());
				_vp.repaint();
			}
		}
	}

	public void UISetNumPeriod() {
		if (_vp.getRNA().get_listeBases().size() != 0) {
			int oldNumPeriod = _vp.getNumPeriod();
			VueNumPeriod vnp = new VueNumPeriod(_vp);
			if (JOptionPane.showConfirmDialog(_vp, vnp.getPanel(),
					"Set new numbering period", JOptionPane.OK_CANCEL_OPTION) != JOptionPane.OK_OPTION) {
				_vp.setNumPeriod(oldNumPeriod);
				_vp.repaint();
			}
		}
	}

	public void UIEditBasePair() {
		if (_vp.isModifiable()) {
			ModeleBase mb = _vp.getRNA().get_listeBases().get(
					_vp.getNearestBase());
			if (mb.getElementStructure() != -1) {
				ModeleStyleBP msbp = mb.getStyleBP();
				ModeleStyleBP.Edge bck5 = msbp.getEdgePartner5();
				ModeleStyleBP.Edge bck3 = msbp.getEdgePartner3();
				ModeleStyleBP.Stericity bcks = msbp.getStericity();

				VueBPType vbpt = new VueBPType(_vp, msbp);

				if (JOptionPane.showConfirmDialog(_vp, vbpt.getPanel(),
						"Set base pair L/W type", JOptionPane.OK_CANCEL_OPTION) != JOptionPane.OK_OPTION) {
					msbp.setEdge5(bck5);
					msbp.setEdge3(bck3);
					msbp.setStericity(bcks);
					_vp.repaint();
				}
			}
		}
	}

	public void UIColorBasePair() {
		if (_vp.isModifiable()) {
			ModeleBase mb = _vp.getRNA().get_listeBases().get(
					_vp.getNearestBase());
			if (mb.getElementStructure() != -1) {
				ModeleStyleBP msbp = mb.getStyleBP();
				Color c = JColorChooser.showDialog(_vp,
						"Choose custom base pair color", msbp.getColor(_vp
								.getConfig()._bondColor));
				if (c != null) {
					msbp.setCustomColor(c);
					_vp.repaint();
				}
			}
		}
	}

	public void UIThicknessBasePair() {
		if (_vp.isModifiable()) {
			ModeleBase mb = _vp.getRNA().get_listeBases().get(
					_vp.getNearestBase());
			if (mb.getElementStructure() != -1) {
				ModeleStyleBP msbp = mb.getStyleBP();
				ArrayList<ModeleStyleBP> bases = new ArrayList<ModeleStyleBP>();
				bases.add(msbp);
				VueBPThickness vbpt = new VueBPThickness(_vp, bases);
				if (JOptionPane.showConfirmDialog(_vp, vbpt.getPanel(),
						"Set base pair(s) thickness",
						JOptionPane.OK_CANCEL_OPTION) != JOptionPane.OK_OPTION) {
					vbpt.restoreThicknesses();
					_vp.repaint();
				}
			}
		}
	}
	

	public void saveToPNG(String filename) throws ExceptionExportFailed {
		VueJPEG jpeg = new VueJPEG(true, false);
		if (JOptionPane.showConfirmDialog(_vp, jpeg.getPanel(),
				"Set resolution", JOptionPane.OK_CANCEL_OPTION) == JOptionPane.OK_OPTION) {
			Double scale = jpeg.getScaleSlider().getValue() / 100.0;
			BufferedImage myImage = new BufferedImage((int) Math.round(_vp
					.getWidth()
					* scale), (int) Math.round(_vp.getHeight() * scale),
					BufferedImage.TRANSLUCENT);
			Graphics2D g2 = myImage.createGraphics();
			AffineTransform AF = new AffineTransform();
			AF.setToScale(scale, scale);
			g2.setTransform(AF);
			_vp.paintComponent(g2,!_vp.getConfig()._drawBackground);
			g2.dispose();
			try {
				ImageIO.write(myImage, "PNG", new File(filename));
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}

	public void saveToJPEG(String filename) throws ExceptionJPEGEncoding,
			ExceptionExportFailed {
		VueJPEG jpeg = new VueJPEG(true, true);
		if (JOptionPane.showConfirmDialog(_vp, jpeg.getPanel(),
				"Set resolution/quality", JOptionPane.OK_CANCEL_OPTION) == JOptionPane.OK_OPTION) {
			Double scale;
			if (jpeg.getScaleSlider().getValue() == 0)
				scale = 1. / 100.;
			else
				scale = jpeg.getScaleSlider().getValue() / 100.;
			BufferedImage myImage = new BufferedImage((int) Math.round(_vp
					.getWidth()
					* scale), (int) Math.round(_vp.getHeight() * scale),
					BufferedImage.TYPE_INT_RGB);
			Graphics2D g2 = myImage.createGraphics();
			AffineTransform AF = new AffineTransform();
			AF.setToScale(scale, scale);
			g2.setTransform(AF);
			_vp.paintComponent(g2);
			try {
				FileImageOutputStream out = new FileImageOutputStream(new File(filename));
				ImageWriter writer = ImageIO.getImageWritersByFormatName("jpeg").next();
				ImageWriteParam params = writer.getDefaultWriteParam();
				params.setCompressionMode(ImageWriteParam.MODE_EXPLICIT);
				params.setCompressionQuality(jpeg.getQualitySlider().getValue() / 100.0f);
				writer.setOutput(out);
				IIOImage myIIOImage = new IIOImage(myImage, null, null);
				writer.write(null, myIIOImage, params);
				out.close();
			} catch (IOException e) {
				throw new ExceptionExportFailed(e.getMessage(), filename);
			}
		}
	}

	public void UIToggleShowNCBP() {
		if (_vp.isModifiable()) {
			_vp.setShowNonCanonicalBP(!_vp.getShowNonCanonicalBP());
			_vp.repaint();
		}
	}

	public void UIToggleColorSpecialBases() {
		_vp.setColorNonStandardBases(!_vp.getColorSpecialBases());
		_vp.repaint();
	}

	public void UIToggleColorGapsBases() {
		_vp.setColorGapsBases(!_vp.getColorGapsBases());
		_vp.repaint();
	}

	public void UIToggleShowNonPlanar() {
		if (_vp.isModifiable()) {
			_vp.setShowNonPlanarBP(!_vp.getShowNonPlanarBP());
			_vp.repaint();
		}
	}

	public void UIToggleShowWarnings() {
		_vp.setShowWarnings(!_vp.getShowWarnings());
		_vp.repaint();
	}

	public void UIPickSpecialBasesColor() {
		Color c = JColorChooser.showDialog(_vp,
				"Choose new special bases color", _vp
						.getNonStandardBasesColor());
		if (c != null) {
			_vp.setNonStandardBasesColor(c);
			_vp.setColorNonStandardBases(true);
			_vp.repaint();
		}
	}

	public void UIPickGapsBasesColor() {
		Color c = JColorChooser.showDialog(_vp, "Choose new gaps bases color",
				_vp.getGapsBasesColor());
		if (c != null) {
			_vp.setGapsBasesColor(c);
			_vp.setColorGapsBases(true);
			_vp.repaint();
		}
	}

	public void UIBaseTypeColor() {
		if (_vp.isModifiable()) {
			new VueBases(_vp, VueBases.KIND_MODE);
		}
	}

	public void UIToggleModifiable() {
		_vp.setModifiable(!_vp.isModifiable());
	}

	public void UIBasePairTypeColor() {
		if (_vp.isModifiable()) {
			new VueBases(_vp, VueBases.COUPLE_MODE);
		}
	}

	public void UIBaseAllColor() {
		if (_vp.isModifiable()) {
			new VueBases(_vp, VueBases.ALL_MODE);
		}
	}

	public void UIAbout() {
		VueAboutPanel about = new VueAboutPanel();
		JOptionPane.showMessageDialog(_vp, about, "About VARNA "
				+ VARNAConfig.MAJOR_VERSION + "." + VARNAConfig.MINOR_VERSION,
				JOptionPane.PLAIN_MESSAGE);
		about.gracefulStop();
	}

	public void UIAutoAnnotateHelices() {
		if (_vp.isModifiable()) {
			_vp.getRNA().autoAnnotateHelices();
			_vp.repaint();
		}
	}

	public void UIAutoAnnotateStrandEnds() {
		if (_vp.isModifiable()) {
			_vp.getRNA().autoAnnotateStrandEnds();
			_vp.repaint();
		}
	}
	
	public void UIAutoAnnotateInteriorLoops() {
		if (_vp.isModifiable()) {
			_vp.getRNA().autoAnnotateInteriorLoops();
			_vp.repaint();
		}
	}

	public void UIAutoAnnotateTerminalLoops() {
		if (_vp.isModifiable()) {
			_vp.getRNA().autoAnnotateTerminalLoops();
			_vp.repaint();
		}
	}

	public void UIAnnotationRemoveFromAnnotation(TextAnnotation textAnnotation) {
		if (_vp.isModifiable()) {
			_vp.set_selectedAnnotation(null);
			_vp.getListeAnnotations().remove(textAnnotation);
			_vp.repaint();
		}
	}

	public void UIAnnotationEditFromAnnotation(TextAnnotation textAnnotation) {
		VueAnnotation vue;
		if (textAnnotation.getType() == TextAnnotation.POSITION)
			vue = new VueAnnotation(_vp, textAnnotation, false);
		else
			vue = new VueAnnotation(_vp, textAnnotation, true, false);
		vue.show();
	}

	public void UIAnnotationAddFromStructure(int type,
			ArrayList<Integer> listeIndex) throws Exception {
		TextAnnotation textAnnot;
		ArrayList<ModeleBase> listeBase;
		VueAnnotation vue;
		switch (type) {
		case TextAnnotation.BASE:
			textAnnot = new TextAnnotation("", _vp.getRNA().get_listeBases()
					.get(listeIndex.get(0)));
			vue = new VueAnnotation(_vp, textAnnot, true);
			vue.show();
			break;
		case TextAnnotation.LOOP:
			listeBase = new ArrayList<ModeleBase>();
			for (Integer i : listeIndex) {
				listeBase.add(_vp.getRNA().get_listeBases().get(i));
			}
			textAnnot = new TextAnnotation("", listeBase, type);
			vue = new VueAnnotation(_vp, textAnnot, true);
			vue.show();
			break;
		case TextAnnotation.HELIX:
			listeBase = new ArrayList<ModeleBase>();
			for (Integer i : listeIndex) {
				listeBase.add(_vp.getRNA().get_listeBases().get(i));
			}
			textAnnot = new TextAnnotation("", listeBase, type);
			vue = new VueAnnotation(_vp, textAnnot, true);
			vue.show();
			break;
		default:
			_vp.errorDialog(new Exception("Unknown structure type"));
			break;
		}
	}

	public void UIAnnotationEditFromStructure(int type,
			ArrayList<Integer> listeIndex) {
		if (_vp.isModifiable()) {
			ModeleBase mb = _vp.getRNA().get_listeBases()
					.get(listeIndex.get(0));
			TextAnnotation ta = _vp.getRNA().getAnnotation(type, mb);
			if (ta != null)
				UIAnnotationEditFromAnnotation(ta);
		}
	}

	public void UIAnnotationRemoveFromStructure(int type,
			ArrayList<Integer> listeIndex) {
		if (_vp.isModifiable()) {
			ModeleBase mb = _vp.getRNA().get_listeBases()
					.get(listeIndex.get(0));
			TextAnnotation ta = _vp.getRNA().getAnnotation(type, mb);
			if (ta != null)
				UIAnnotationRemoveFromAnnotation(ta);
		}
	}

	public void UIAnnotationsAddPosition(int x, int y) {
		if (_vp.isModifiable()) {
			Point2D.Double p = _vp.panelToLogicPoint(new Point2D.Double(x, y));
			VueAnnotation annotationAdd = new VueAnnotation(_vp, (int) p.x,
					(int) p.y);
			annotationAdd.show();
		}
	}

	public void UIAnnotationsAddBase(int x, int y) {
		if (_vp.isModifiable()) {
			ModeleBase mb = _vp.getBaseAt(new Point2D.Double(x, y));
			if(mb!=null)
			{
				_vp.setSelectedBase(mb.getIndex());
				TextAnnotation textAnnot = new TextAnnotation("", mb);
				VueAnnotation annotationAdd = new VueAnnotation(_vp,  textAnnot,true);
				annotationAdd.show();
			}
		}
	}

	public void UIAnnotationsAddLoop(int x, int y) {
		if (_vp.isModifiable()) {
			try {
			ModeleBase mb = _vp.getBaseAt(new Point2D.Double(x, y));
			if(mb!=null)
			{
				Vector<Integer> v = _vp.getRNA().getLoopBases(mb.getIndex());
				ArrayList<ModeleBase> mbs = _vp.getRNA().getBasesAt(v); 
				TextAnnotation textAnnot;
					textAnnot = new TextAnnotation("", mbs,TextAnnotation.LOOP);
				_vp.setSelection(mbs);
				VueAnnotation annotationAdd = new VueAnnotation(_vp,  textAnnot,true);
				annotationAdd.show();
			}
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}
	
	private ArrayList<ModeleBase> extractMaxContiguousPortion(ArrayList<ModeleBase> m)
	{
		ModeleBase[] tab = new ModeleBase[_vp.getRNA().getSize()];
		for(int i=0;i<tab.length;i++)
		{ tab[i] = null; }
		for(ModeleBase mb : m)
		{ tab[mb.getIndex()] = mb; }
		ArrayList<ModeleBase> best = new ArrayList<ModeleBase>();
		ArrayList<ModeleBase> current = new ArrayList<ModeleBase>();
		for(int i=0;i<tab.length;i++)
		{ 
			if (tab[i] != null)
			{  current.add(tab[i]);  }
			else
			{    
				if (current.size()>best.size())
				  best = current;
				current = new ArrayList<ModeleBase>();
			}
		}
		if (current.size()>best.size())
		{
			best = current;
		}
		return best;
	}


	
	public void UIAnnotationsAddRegion(int x, int y) {
		if (_vp.isModifiable()) {
			ArrayList<ModeleBase> mb = _vp.getSelection().getBases();
			if (mb.size()==0)
			{
				ModeleBase m = _vp.getBaseAt(new Point2D.Double(x, y));
				mb.add(m);
			}
			mb = extractMaxContiguousPortion(extractMaxContiguousPortion(mb));
			_vp.setSelection(mb);
				HighlightRegionAnnotation regionAnnot = new HighlightRegionAnnotation(mb);
				_vp.addHighlightRegion(regionAnnot);
				VueHighlightRegionEdit annotationAdd = new VueHighlightRegionEdit(_vp,regionAnnot);
				if (!annotationAdd.show())
				{  _vp.removeHighlightRegion(regionAnnot);  }
				_vp.clearSelection();
		}
	}

	public void UIAnnotationsAddChemProb(int x, int y) {
		if (_vp.isModifiable() && _vp.getRNA().getSize()>1) {
			Point2D.Double p = _vp.panelToLogicPoint(new Point2D.Double(x, y));
			ModeleBase m1 = _vp.getBaseAt(new Point2D.Double(x, y));
			ModeleBase best = null;
			if (m1.getIndex()-1>=0)
			{ best = _vp.getRNA().getBaseAt(m1.getIndex()-1);}
			if (m1.getIndex()+1<_vp.getRNA().getSize())
			{
				ModeleBase m2 = _vp.getRNA().getBaseAt(m1.getIndex()+1);
				if (best==null)
				{
					best = m2;
				}
				else
				{
					if (best.getCoords().distance(p)>m2.getCoords().distance(p))
					{
						best = m2;
					}
				}
			}
			ArrayList<ModeleBase> tab = new ArrayList<ModeleBase>();
			tab.add(m1);
			tab.add(best);
			_vp.setSelection(tab);
				ChemProbAnnotation regionAnnot = new ChemProbAnnotation(m1,best);
				_vp.getRNA().addChemProbAnnotation(regionAnnot);
				VueChemProbAnnotation annotationAdd = new VueChemProbAnnotation(_vp,regionAnnot);
				if (!annotationAdd.show())
				{  _vp.getRNA().removeChemProbAnnotation(regionAnnot);  }
				_vp.clearSelection();
		}
	}

	
	public void UIAnnotationsAddHelix(int x, int y) {
		if (_vp.isModifiable()) {
			try {
			ModeleBase mb = _vp.getBaseAt(new Point2D.Double(x, y));
			if(mb!=null)
			{
				ArrayList<Integer> v = _vp.getRNA().findHelix(mb.getIndex());
				ArrayList<ModeleBase> mbs = _vp.getRNA().getBasesAt(v); 
				TextAnnotation textAnnot;
					textAnnot = new TextAnnotation("", mbs,TextAnnotation.HELIX);
				_vp.setSelection(mbs);
				VueAnnotation annotationAdd = new VueAnnotation(_vp,  textAnnot,true);
				annotationAdd.show();
			}
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}

	public void UIToggleGaspinMode() {
		if (_vp.isModifiable()) {
			_vp.toggleDrawOutlineBase();
			_vp.toggleFillBase();
			_vp.repaint();
		}
	}

	
	public void UIAnnotationsAdd() {
		if (_vp.isModifiable()) {
			VueAnnotation annotationAdd = new VueAnnotation(_vp);
			annotationAdd.show();
		}
	}

	public void UIAnnotationsRemove() {
		if (_vp.isModifiable()) {
			new VueListeAnnotations(_vp, VueListeAnnotations.REMOVE);
		}
	}

	public void UIAnnotationsEdit() {
		if (_vp.isModifiable()) {
			new VueListeAnnotations(_vp, VueListeAnnotations.EDIT);
		}
	}
}
