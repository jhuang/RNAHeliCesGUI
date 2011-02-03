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
package fr.orsay.lri.varna.controlers;

import java.awt.Component;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Vector;

import javax.swing.JMenu;
import javax.swing.JMenuItem;
import javax.swing.event.PopupMenuEvent;
import javax.swing.event.PopupMenuListener;

import fr.orsay.lri.varna.VARNAPanel;
import fr.orsay.lri.varna.models.annotations.TextAnnotation;
import fr.orsay.lri.varna.models.rna.ModeleBase;
import fr.orsay.lri.varna.models.rna.ModeleBaseNucleotide;
import fr.orsay.lri.varna.models.rna.ModeleBasesComparison;
import fr.orsay.lri.varna.models.rna.RNA;


/**
 * Controller of the mouse click
 * 
 * @author darty
 * 
 */
public class ControleurClicMovement implements MouseListener,
		MouseMotionListener, PopupMenuListener {
	private VARNAPanel _vp;
	private boolean _presenceMenuSelection;
	private JMenu _submenuSelection;
	public Point _spawnPoint;	
	public Point _selectionInitialPoint;
	public Point _selectionCurrentPoint;

	public static final double MIN_SELECTION_DISTANCE = 40.0;
	public static final double HYSTERESIS_DISTANCE = 10.0;
	
	
	public enum MouseStates {
		NONE,
		MOVE_ELEMENT,
		SELECT_ELEMENT,
		SELECT_REGION_OR_UNSELECT,
		SELECT_REGION,
		POPUP_MENU,
		MOVE_ANNOTATION,
	};
	private MouseStates _currentState = MouseStates.NONE;
	


	public ControleurClicMovement(VARNAPanel _vuep) {
		_vp = _vuep;
		_vp.getPopup().addPopupMenuListener(this);
		_presenceMenuSelection = false;
	}

	public void mouseClicked(MouseEvent arg0) {
	}

	public void mouseEntered(MouseEvent arg0) {
	}

	public void mouseExited(MouseEvent arg0) {
	}

	public void mousePressed(MouseEvent arg0) 
	{
		boolean button1 = (arg0.getButton() == MouseEvent.BUTTON1);
		boolean button2 = (arg0.getButton() == MouseEvent.BUTTON2);
		boolean button3 = (arg0.getButton() == MouseEvent.BUTTON3);
		boolean shift = arg0.isShiftDown();
		boolean ctrl = arg0.isControlDown();
		boolean alt = arg0.isAltDown();
		if (button1 && !ctrl && !alt && !shift)
		{
			if (_vp.isModifiable()) 
			{
				_currentState = MouseStates.MOVE_ELEMENT;
				if (_vp.getRealCoords() != null
						&& _vp.getRealCoords().length != 0
						&& _vp.getRNA().get_listeBases().size() != 0) 
				{
					int selectedIndex = this.getNearestBaseIndex(arg0,false,_vp.getRNA().get_drawMode()==RNA.DRAW_MODE_RADIATE);
					TextAnnotation selectedAnnotation = this.getNearestAnnotation(arg0);
					if (selectedIndex !=-1)
					{ 
						_vp.setSelectedBase(selectedIndex);
						if (_vp.getRNA().get_drawMode() == RNA.DRAW_MODE_RADIATE) {
							_vp.highlightSelectedStem();
						} else {
							_vp.highlightSelectedBase();
						}
					}
					else
					{
						if (selectedAnnotation != null)
						{
							_currentState = MouseStates.MOVE_ANNOTATION;
							_vp.set_selectedAnnotation(selectedAnnotation);
							_vp.highlightSelectedAnnotation();
						}
						else
						{
							_currentState = MouseStates.SELECT_REGION_OR_UNSELECT;
							_selectionInitialPoint = new Point(arg0.getX(),arg0.getY());
							_selectionCurrentPoint = new Point(_selectionInitialPoint);
						}
					}
				}
			}
		}
		if (button1 && !ctrl && !alt && shift)
		{
			_currentState = MouseStates.SELECT_ELEMENT;
			_selectionInitialPoint = new Point(arg0.getX(),arg0.getY());
			_selectionCurrentPoint = new Point(_selectionInitialPoint);
		}
		if (button1 && ctrl && !alt && !shift)
		{
			_currentState = MouseStates.SELECT_ELEMENT;
			_selectionInitialPoint = new Point(arg0.getX(),arg0.getY());
			_selectionCurrentPoint = new Point(_selectionInitialPoint);
		}
		if (button3) 
		{
			_currentState = MouseStates.POPUP_MENU;
			// si le menu sur la selection est deja presente dans le menu popup
			// on la retire ainsi que le separateur
			if (_presenceMenuSelection) {
				_vp.getPopupMenu().removeSelectionMenu();
			}

			// si il y a des bases
			if (_vp.getRealCoords() != null
					&& _vp.getRNA().get_listeBases().size() != 0) {
				// on calcul la base la plus proche
				updateNearestBase(arg0);
				// on insere dans le menu les nouvelles options
				addMenu(arg0);
				if (_vp.get_selectedAnnotation() != null)
					_vp.highlightSelectedAnnotation();
			}
			// affichage du popup menu
			if (_vp.getRNA().get_drawMode() == RNA.DRAW_MODE_LINEAR) {
				_vp.getPopup().get_rotation().setEnabled(false);
			} else {
				_vp.getPopup().get_rotation().setEnabled(true);
			}
			_vp.getPopup().updateDialog();
			_vp.getPopup().show(_vp, arg0.getX(), arg0.getY());
		}
		_vp.repaint();
	}

	public void mouseReleased(MouseEvent arg0) {
		if (arg0.getButton() == MouseEvent.BUTTON1)
		{
			if (_currentState == MouseStates.MOVE_ELEMENT)
			{
				_vp.clearSelection();
				_vp.setSelectedBase(-1);
				_vp.removeSelectedAnnotation();
			}
			else if (_currentState == MouseStates.SELECT_REGION_OR_UNSELECT)
			{
				_vp.clearSelection();
				_vp.setSelectedBase(-1);
				_vp.removeSelectedAnnotation();
			}
			else if (_currentState == MouseStates.SELECT_ELEMENT)
			{
				if (_vp.getRealCoords() != null
						&& _vp.getRealCoords().length != 0
						&& _vp.getRNA().get_listeBases().size() != 0) 
					{
						int selectedIndex = this.getNearestBaseIndex(arg0,false,false);
						if (selectedIndex !=-1)
						{ 
							_vp.toggleSelection(selectedIndex);
						}
					}
				_vp.setSelectedBase(-1);
				_vp.removeSelectedAnnotation();
			}
			else if (_currentState == MouseStates.SELECT_REGION)
			{
			  _vp.removeSelectionRectangle();
			}

		}
		_vp.repaint();				
		_currentState = MouseStates.NONE;
	}

	public void mouseDragged(MouseEvent me) {
		if (_currentState == MouseStates.MOVE_ELEMENT)
		{
				// si on deplace la souris et qu'une base est selectionnÃ©e
				if (_vp.getSelectedBaseIndex() != -1) {
					if (_vp.getRNA().get_drawMode() == RNA.DRAW_MODE_CIRCULAR
							|| _vp.getRNA().get_drawMode() == RNA.DRAW_MODE_LINEAR
							|| _vp.getRNA().get_drawMode() == RNA.DRAW_MODE_NAVIEW
							|| _vp.getRNA().get_drawMode() == RNA.DRAW_MODE_TEMPLATE) {
						// dans le cas circulaire naview ou line on deplace une base
						moveSingleAtom(_vp.getSelectedBaseIndex(), me.getX(), me.getY());
					} else if (_vp.getRNA().get_drawMode() == RNA.DRAW_MODE_RADIATE) {
						// dans le cas radiale on deplace une helice
						moveHelixAtom(_vp.getSelectedBaseIndex(), me.getX(), me.getY());
					}
					_vp.repaint();
				}
		}
		else if (_currentState == MouseStates.MOVE_ANNOTATION)
		{
			if (_vp.get_selectedAnnotation()!=null)
			{
				Point2D.Double p = _vp.panelToLogicPoint(new Point2D.Double(me.getX(), me.getY()));
				_vp.get_selectedAnnotation().setAncrage(p.x,p.y);
				_vp.repaint();
			}			
		}
		else if ((_currentState == MouseStates.SELECT_ELEMENT)||(_currentState == MouseStates.SELECT_REGION_OR_UNSELECT))
		{
			if (_selectionInitialPoint.distance(me.getX(),me.getY())>HYSTERESIS_DISTANCE)
				_currentState = MouseStates.SELECT_REGION;
		}
		else if (_currentState == MouseStates.SELECT_REGION)
		{
			_selectionCurrentPoint = new Point(me.getX(),me.getY());
			int minx = Math.min(_selectionCurrentPoint.x, _selectionInitialPoint.x);
			int miny = Math.min(_selectionCurrentPoint.y, _selectionInitialPoint.y);
			int maxx = Math.max(_selectionCurrentPoint.x, _selectionInitialPoint.x);
			int maxy = Math.max(_selectionCurrentPoint.y, _selectionInitialPoint.y);
			_vp.setSelectionRectangle(new Rectangle(minx,miny,maxx-minx,maxy-miny));
		}
		
	}

	private void addMenu(MouseEvent arg0) {
		// creation du menu
		_submenuSelection = new JMenu("Selection");
		addCurrent();
		// ajout des option sur base
		addMenuBase();
		// ajout des option sur paire de base
		if (_vp.getRNA().get_listeBases().get(_vp.getNearestBase())
				.getElementStructure() != -1) {
			addMenuBasePair();
		}

		// detection renflement
		detectBulge();
		// detection 3'
		detect3Prime();
		// detection 5'
		detect5Prime();
		// detection boucle
		detectLoop();
		// detection d'helice
		detectHelix();
		// detection tige
		detectStem();
		// Ajout de toutes bases
		addAllBase();
		// detection d'annotation
		detectAnnotation(arg0);

		_vp.getPopup().addSelectionMenu(_submenuSelection);
		_presenceMenuSelection = true;
	}

	private void detectAnnotation(MouseEvent arg0) {
		if (_vp.getListeAnnotations().size() != 0) {
			double dist = Double.MAX_VALUE;
			double d2;
			Point2D.Double position;
			for (TextAnnotation textAnnot : _vp.getListeAnnotations()) {
				// calcul de la distance
				position = textAnnot.getCenterPosition();
				position = _vp.transformCoord(position);
				d2 = Math.sqrt(Math.pow((position.x - arg0.getX()), 2)
						+ Math.pow((position.y - arg0.getY()), 2));
				// si la valeur est inferieur au minimum actuel
				if (dist > d2) {
					_vp.set_selectedAnnotation(textAnnot);
					dist = d2;
				}
			}
			_submenuSelection.addSeparator();
			_vp.getPopup().addAnnotationMenu(_submenuSelection,true);
		}
	}

	private void detectBulge() {
		int indiceB = _vp.getNearestBase();
		ArrayList<Integer> indices = _vp.getRNA().findBulge(indiceB);
		if ((indices.size() > 0)
				&& (_vp.getRNA().getHelixCountOnLoop(_vp.getNearestBase()) == 2)) {
			JMenu submenuBulge = new JMenu("Bulge");
			submenuBulge.addChangeListener(new ControleurSelectionHighlight(
					new Vector<Integer>(indices), _vp, submenuBulge));
			submenuBulge.setActionCommand("bulge");
			if (!_vp.isModifiable())
				submenuBulge.setEnabled(false);
			_vp.getPopupMenu().addColorOptions(submenuBulge);
			_submenuSelection.add(submenuBulge);
		}
	}

	private void detectHelix() {
		int indiceH = _vp.getNearestBase();
		ArrayList<Integer> indices = _vp.getRNA().findHelix(indiceH);
		if (indices.size() != 0) {
			// ajout menu helice
			JMenu submenuHelix = new JMenu("Helix");
			submenuHelix.addChangeListener(new ControleurSelectionHighlight(
					new Vector<Integer>(indices), _vp, submenuHelix));
			submenuHelix.setActionCommand("helix");
			if (!_vp.isModifiable())
				submenuHelix.setEnabled(false);
			_vp.getPopupMenu().addColorOptions(submenuHelix);
			submenuHelix.addSeparator();
			_vp.getPopupMenu().addAnnotationMenu(submenuHelix);
			_submenuSelection.add(submenuHelix);
		}
	}

	private void detectStem() {
		int indiceS = _vp.getNearestBase();
		ArrayList<Integer> indices = _vp.getRNA().findStem(indiceS);
		if (indices.size() > 0) {
			JMenu submenuStem = new JMenu("Stem");
			submenuStem.addChangeListener(new ControleurSelectionHighlight(
					new Vector<Integer>(indices), _vp, submenuStem));
			submenuStem.setActionCommand("stem");
			if (!_vp.isModifiable())
				submenuStem.setEnabled(false);
			_vp.getPopupMenu().addColorOptions(submenuStem);
			_submenuSelection.add(submenuStem);
		}
	}

	private void detect3Prime() {
		// detection 3'
		int indice3 = _vp.getNearestBase();
		ArrayList<Integer> indices = _vp.getRNA().find3Prime(indice3);
		if (indices.size() != 0) {
			JMenu submenu3Prime = new JMenu("3'");
			submenu3Prime.addChangeListener(new ControleurSelectionHighlight(
					new Vector<Integer>(indices), _vp, submenu3Prime));
			submenu3Prime.setActionCommand("3'");
			if (!_vp.isModifiable())
				submenu3Prime.setEnabled(false);
			_vp.getPopupMenu().addColorOptions(submenu3Prime);
			_submenuSelection.add(submenu3Prime);
		}
	}

	private void detect5Prime() {
		int indice5 = _vp.getNearestBase();
		ArrayList<Integer> indices = _vp.getRNA().find5Prime(indice5);
		if (indices.size() != 0) {
			JMenu submenu5Prime = new JMenu("5'");
			submenu5Prime.addChangeListener(new ControleurSelectionHighlight(
					new Vector<Integer>(indices), _vp, submenu5Prime));
			submenu5Prime.setActionCommand("5'");
			if (!_vp.isModifiable())
				submenu5Prime.setEnabled(false);
			_vp.getPopupMenu().addColorOptions(submenu5Prime);
			_submenuSelection.add(submenu5Prime);
		}
	}

	private void detectLoop() {
		int indexL = _vp.getNearestBase();
		if (_vp.getRNA().get_listeBases().get(indexL).getElementStructure() == -1) {
			ArrayList<Integer> listLoop = _vp.getRNA().findLoop(indexL);
			JMenu submenuLoop = new JMenu("Loop");
			submenuLoop.addChangeListener(new ControleurSelectionHighlight(
					listLoop, _vp, submenuLoop));
			submenuLoop.setActionCommand("loop1");
			if (!_vp.isModifiable())
				submenuLoop.setEnabled(false);
			_vp.getPopupMenu().addColorOptions(submenuLoop);
			submenuLoop.addSeparator();
			_vp.getPopupMenu().addAnnotationMenu(submenuLoop);
			_submenuSelection.add(submenuLoop);
		} else {
			ArrayList<Integer> listLoop1 = _vp.getRNA().findLoopForward(indexL);
			if (listLoop1.size() > 0) {
				JMenu submenuLoop1 = new JMenu("Forward loop");
				submenuLoop1
						.addChangeListener(new ControleurSelectionHighlight(
								listLoop1, _vp, submenuLoop1));
				submenuLoop1.setActionCommand("loop1");
				if (!_vp.isModifiable())
					submenuLoop1.setEnabled(false);
				_vp.getPopupMenu().addColorOptions(submenuLoop1);
				submenuLoop1.addSeparator();
				_vp.getPopupMenu().addAnnotationMenu(submenuLoop1);
				_submenuSelection.add(submenuLoop1);
			}
			ArrayList<Integer> listLoop2 = _vp.getRNA()
					.findLoopBackward(indexL);
			if (listLoop2.size() > 0) {
				JMenu submenuLoop2 = new JMenu("Backward loop");
				submenuLoop2
						.addChangeListener(new ControleurSelectionHighlight(
								listLoop2, _vp, submenuLoop2));
				submenuLoop2.setActionCommand("loop2");
				if (!_vp.isModifiable())
					submenuLoop2.setEnabled(false);
				_vp.getPopupMenu().addColorOptions(submenuLoop2);
				submenuLoop2.addSeparator();
				_vp.getPopupMenu().addAnnotationMenu(submenuLoop2);
				_submenuSelection.add(submenuLoop2);
			}
		}
	}

	private void addCurrent() {
		Collection<? extends ModeleBase> mbs = _vp.getSelection().getBases();
		if (mbs.size()>0)
		{
		JMenu submenuAll = new JMenu("Current");
		submenuAll.addChangeListener(new ControleurSelectionHighlight(
				mbs, _vp, submenuAll));
		submenuAll.setActionCommand("current");
		if (!_vp.isModifiable())
			submenuAll.setEnabled(false);
		_vp.getPopupMenu().addColorOptions(submenuAll);
		_submenuSelection.add(submenuAll);
		}
	}
	
	
	private void addMenuBase() {
		JMenu submenuBase = new JMenu();
		ModeleBase mb = _vp.getRNA().get_listeBases().get(_vp.getNearestBase());
		if (_vp.getRNA().is_comparisonMode()) {
			submenuBase.setText("Base #" + (mb.getBaseNumber()) + ":"
					+ ((ModeleBasesComparison) mb).getBases());
		} else {
			submenuBase.setText("Base #" + (mb.getBaseNumber()) + ":"
					+ ((ModeleBaseNucleotide) mb).get_c());
		}
		submenuBase.addChangeListener(new ControleurSelectionHighlight(mb
				.getIndex(), _vp, submenuBase));
		submenuBase.setActionCommand("base");
		// option disponible seulement en mode modifiable
		if (!_vp.isModifiable())
			submenuBase.setEnabled(false);

		JMenuItem baseChar = new JMenuItem("Edit base");
		baseChar.setActionCommand("baseChar");
		baseChar.addActionListener(_vp.getPopupMenu().get_controleurMenu());
		submenuBase.add(baseChar);
		_vp.getPopupMenu().addColorOptions(submenuBase);
		submenuBase.addSeparator();
		_vp.getPopupMenu().addAnnotationMenu(submenuBase);
		_submenuSelection.add(submenuBase);
	}

	private void addAllBase() {
		ArrayList<Integer> indices = _vp.getRNA().findAll();
		JMenu submenuAll = new JMenu("All");
		submenuAll.addChangeListener(new ControleurSelectionHighlight(
				new Vector<Integer>(indices), _vp, submenuAll));
		submenuAll.setActionCommand("all");
		if (!_vp.isModifiable())
			submenuAll.setEnabled(false);
		_vp.getPopupMenu().addColorOptions(submenuAll);
		_submenuSelection.add(submenuAll);
	}

	private void addMenuBasePair() {
		if (!_vp.getRNA().is_comparisonMode()) {
			int indiceBP = _vp.getNearestBase();
			ArrayList<Integer> indices = _vp.getRNA().findPair(indiceBP);
			ModeleBaseNucleotide base = ((ModeleBaseNucleotide) _vp.getRNA()
					.get_listeBases().get(_vp.getNearestBase()));
			if (base.getElementStructure() != -1) {
				JMenu submenuBasePair = new JMenu();
				ModeleBaseNucleotide partner = ((ModeleBaseNucleotide) _vp
						.getRNA().get_listeBases().get(
								base.getElementStructure()));
				submenuBasePair
						.addChangeListener(new ControleurSelectionHighlight(
								indices, _vp, submenuBasePair));
				submenuBasePair.setText("Base pair #("
						+ (Math.min(base.getBaseNumber(), partner
								.getBaseNumber()))
						+ ","
						+ (Math.max(base.getBaseNumber(), partner
								.getBaseNumber())) + ")");
				submenuBasePair.setActionCommand("bp");
				// option disponible seulement en mode modifiable
				if (!_vp.isModifiable())
					submenuBasePair.setEnabled(false);

				JMenuItem basepair = new JMenuItem("Edit BP");
				basepair.setActionCommand("basepair");
				basepair.addActionListener(_vp.getPopupMenu()
						.get_controleurMenu());

				_vp.getPopupMenu().addColorOptions(submenuBasePair);
				Component[] comps = submenuBasePair.getMenuComponents();
				int offset = -1;
				for (int i = 0; i < comps.length; i++) {
					Component c = comps[i];
					if (c instanceof JMenuItem) {
						JMenuItem jmi = (JMenuItem) c;
						if (jmi.getActionCommand().contains(",BPColor")) {
							offset = i;
						}
					}
				}
				if (offset != -1) {
					submenuBasePair.insert(basepair, offset);
				} else {
					submenuBasePair.add(basepair);
				}
				_submenuSelection.add(submenuBasePair);
			}
		}
	}
	private ModeleBase getNearestBase(MouseEvent arg0, boolean always, boolean onlyPaired) {
		int i = getNearestBaseIndex(arg0, always,onlyPaired);
		if (i==-1) return null;
    	return  _vp.getRNA().get_listeBases().get(i);
	}

	private int getNearestBaseIndex(MouseEvent arg0) {
		return getNearestBaseIndex(arg0, false ,false);
	}
	private int getNearestBaseIndex(MouseEvent arg0, boolean always, boolean onlyPaired) {
		double d2, dist = Double.MAX_VALUE;
		int mb = -1;
		for (int i = 0; i < _vp.getRealCoords().length; i++) {
			if (!onlyPaired || (_vp.getRNA().get_listeBases().get(i).getElementStructure()!=-1))
			{d2 = Math.sqrt(Math
					.pow((_vp.getRealCoords()[i].x - arg0.getX()), 2)
					+ Math.pow((_vp.getRealCoords()[i].y - arg0.getY()), 2));
			if ((dist > d2) && ((d2<_vp.getScaleFactor()*MIN_SELECTION_DISTANCE) || always)) {
				dist = d2;
				mb = i;
			}
			}
		}
		return mb;
	}	

	private void updateNearestBase(MouseEvent arg0) {
		int i = getNearestBaseIndex(arg0,true,false);
		if (i!=-1)
			_vp.setNearestBase(i);
	}


	private TextAnnotation getNearestAnnotation(MouseEvent arg0) {
		TextAnnotation t = null;
		if (_vp.getListeAnnotations().size() != 0) {
			double dist = Double.MAX_VALUE;
			double d2;
			Point2D.Double position;
			for (TextAnnotation textAnnot : _vp.getListeAnnotations()) {
				// calcul de la distance
				position = textAnnot.getCenterPosition();
				position = _vp.transformCoord(position);
				d2 = Math.sqrt(Math.pow((position.x - arg0.getX()), 2)
						+ Math.pow((position.y - arg0.getY()), 2));
				// si la valeur est inferieur au minimum actuel
				if ((dist > d2)&& (d2<_vp.getScaleFactor()*MIN_SELECTION_DISTANCE)) {
					t = textAnnot;
					dist = d2;
				}
			}
		}
		return t;
	}	

	
	public void mouseMoved(MouseEvent arg0) {
		_vp.setLastSelectedPosition(new Point2D.Double(arg0.getX(),arg0.getY()));
	}

	/**
	 * Move a base of the rna
	 * 
	 * @param index
	 *            :the index of the base to move in the base list
	 * @param x
	 *            :the new x coordinate
	 * @param y
	 *            :the new y coordinate
	 */
	private void moveSingleAtom(int index, int x, int y) {
		if (_vp.isModifiable() && (index >= 0)
				&& (index < _vp.getRNA().get_listeBases().size())) {
			Point2D.Double pr = _vp.getRealCoords()[index];
			Point2D.Double pold = _vp.getRNA().getCoords(index);
			double dx = ((x - pr.x) / _vp.getScaleFactor());
			double dy = ((y - pr.y) / _vp.getScaleFactor());

			_vp.getRNA().setCoord(index, pold.x + dx, pold.y + dy);
			if (_vp.getRNA().get_drawMode() == RNA.DRAW_MODE_LINEAR
					&& _vp.getRNA().get_listeBases().get(index)
							.getElementStructure() != -1) {
				index = _vp.getRNA().get_listeBases().get(index)
						.getElementStructure();
				pold = _vp.getRNA().getCoords(index);
				_vp.getRNA().setCoord(index, pold.x + dx, pold.y + dy);
			}
		}
	}

	private Point2D.Double project(Point2D.Double O, Point2D.Double Ox,
			Point2D.Double C) {
		Point2D.Double OC = new Point2D.Double(C.x - O.x, C.y - O.y);
		// Projection of OC on OI => OX
		double normOX = (Ox.x * OC.x + Ox.y * OC.y);
		Point2D.Double OX = new Point2D.Double((normOX * Ox.x), (normOX * Ox.y));
		// Portion of OC orthogonal to Ox => XC
		Point2D.Double XC = new Point2D.Double(OC.x - OX.x, OC.y - OX.y);
		// Reflexive image of C with respect to Ox => CP
		Point2D.Double OCP = new Point2D.Double(OX.x - XC.x, OX.y - XC.y);
		Point2D.Double CP = new Point2D.Double(O.x + OCP.x, O.y + OCP.y);
		return CP;
	}

	/**
	 * Flip an helix around its supporting base
	 */
	private void flipHelix(int hBeg, int hEnd) {
		Point2D.Double A = _vp.getRNA().getCoords(hBeg);
		Point2D.Double B = _vp.getRNA().getCoords(hEnd);
		Point2D.Double AB = new Point2D.Double(B.x - A.x, B.y - A.y);
		double normAB = Math.sqrt(AB.x * AB.x + AB.y * AB.y);
		// Creating a coordinate system centered on A and having
		// unit x-vector Ox.
		Point2D.Double O = A;
		Point2D.Double Ox = new Point2D.Double(AB.x / normAB, AB.y / normAB);
		for (int i = hBeg + 1; i < hEnd; i++) {
			Point2D.Double P = _vp.getRNA().getCoords(i);
			_vp.getRNA().setCoord(i, project(O, Ox, P));
			Point2D.Double Center = _vp.getRNA().getCenter(i);
			_vp.getRNA().setCenter(i, project(O, Ox, Center));
		}
	}

	/**
	 * Tests if an helix needs to be flipped.
	 */
	boolean shouldFlip(int index, int x, int y) {
		Point h = _vp.getRNA().getHelixInterval(index);

		Point2D.Double P = _vp.panelToLogicPoint(new Point2D.Double(x, y));
		Point2D.Double A = _vp.getRNA().getCoords(h.x);
		Point2D.Double B = _vp.getRNA().getCoords(h.y);
		Point2D.Double C = _vp.getRNA().getCoords(h.x + 1);
		// Creating a vector that is orthogonal to AB
		Point2D.Double hAB = new Point2D.Double(B.y - A.y, -(B.x - A.x));
		Point2D.Double AC = new Point2D.Double(C.x - A.x, C.y - A.y);
		Point2D.Double AP = new Point2D.Double(P.x - A.x, P.y - A.y);
		double signC = (hAB.x * AC.x + hAB.y * AC.y);
		double signP = (hAB.x * AP.x + hAB.y * AP.y);
		// Now, the product signC*signP is negative iff the mouse and the first
		// base inside
		// the helix are on different sides of the end of the helix => Flip the
		// helix!
		return (signC * signP < 0.0);
	}

	private boolean testDirectionality(int i, int j, int k) {

		// Which direction are we heading toward?
		Point2D.Double pi = _vp.getRNA().getCoords(i);
		Point2D.Double pj = _vp.getRNA().getCoords(j);
		Point2D.Double pk = _vp.getRNA().getCoords(k);
		return testDirectionality(pi, pj, pk);
	}

	public static boolean testDirectionality(Point2D.Double pi, Point2D.Double pj, Point2D.Double pk) {

		// Which direction are we heading toward?
		double test = (pj.x - pi.x) * (pk.y - pj.y) - (pj.y - pi.y)
				* (pk.x - pj.x);
		return test < 0.0;
	}

	/**
	 * Move a helix of the rna
	 * 
	 * @param index
	 *            :the index of the selected base
	 * @param x
	 *            :the new x coordinate
	 * @param y
	 *            :the new y coordinate
	 */
	private void moveHelixAtom(int index, int x, int y) {
		if (_vp.isModifiable() && (index >= 0)
				&& (index < _vp.getRNA().get_listeBases().size())) {
			int indexTo = _vp.getRNA().get_listeBases().get(index)
					.getElementStructure();
			Point h = _vp.getRNA().getHelixInterval(index);
			Point ml = _vp.getRNA().getMultiLoop(h.x);
			int i = ml.x;
			if (indexTo != -1) {
				if (i == 0) {
					if (shouldFlip(index, x, y)) {
						flipHelix(h.x, h.y);
					}
				} else {

					int prevIndex = h.x;
					int nextIndex = h.y;
					while (i <= ml.y) {
						int j = _vp.getRNA().get_listeBases().get(i)
								.getElementStructure();
						if ((j != -1) && (i < h.x)) {
							prevIndex = i;
						}
						if ((j != -1) && (i > h.y) && (nextIndex == h.y)) {
							nextIndex = i;
						}
						if ((j > i) && (j < ml.y)) {
							i = _vp.getRNA().get_listeBases().get(i)
									.getElementStructure();
						} else {
							i++;
						}
					}
					Point2D.Double p = new Point2D.Double(x, y);
					Point2D.Double oldPos = _vp.getRNA().getCoords(index);
					Point2D.Double limitLoopLeft, limitLoopRight, limitLeft, limitRight, helixStart, helixStop;
					boolean isDirect = testDirectionality(ml.x, ml.y, h.x);
					if (isDirect) {
						limitLoopLeft = _vp.getRNA().getCoords(ml.y);
						limitLoopRight = _vp.getRNA().getCoords(ml.x);
						limitLeft = _vp.getRNA().getCoords(prevIndex);
						limitRight = _vp.getRNA().getCoords(nextIndex);
						helixStart = _vp.getRNA().getCoords(h.x);
						helixStop = _vp.getRNA().getCoords(h.y);
					} else {
						limitLoopLeft = _vp.getRNA().getCoords(ml.x);
						limitLoopRight = _vp.getRNA().getCoords(ml.y);
						limitLeft = _vp.getRNA().getCoords(nextIndex);
						limitRight = _vp.getRNA().getCoords(prevIndex);
						helixStart = _vp.getRNA().getCoords(h.y);
						helixStop = _vp.getRNA().getCoords(h.x);
					}

					Point2D.Double center = _vp.getRNA().get_listeBases().get(
							h.x).getCenter();
					Point2D.Double newPos = _vp.panelToLogicPoint(p);
					double base = (computeAngle(center, limitLoopRight) + computeAngle(
							center, limitLoopLeft)) / 2.0;
					double pLimR = computeAngle(center, limitLeft) - base;
					double pHelR = computeAngle(center, helixStart) - base;
					double pNew = computeAngle(center, newPos) - base;
					double pOld = computeAngle(center, oldPos) - base;
					double pHelL = computeAngle(center, helixStop) - base;
					double pLimL = computeAngle(center, limitRight) - base;

					while (pLimR < 0.0)
						pLimR += 2.0 * Math.PI;
					while (pHelR < pLimR)
						pHelR += 2.0 * Math.PI;
					while ((pNew < pHelR))
						pNew += 2.0 * Math.PI;
					while ((pOld < pHelR))
						pOld += 2.0 * Math.PI;
					while ((pHelL < pOld))
						pHelL += 2.0 * Math.PI;
					while ((pLimL < pHelL))
						pLimL += 2.0 * Math.PI;

					double minDelta = normalizeAngle((pLimR - pHelR) + 0.25);
					double maxDelta = normalizeAngle((pLimL - pHelL) - 0.25);
					while (maxDelta < minDelta)
						maxDelta += 2.0 * Math.PI;
					double delta = normalizeAngle(pNew - pOld);
					while (delta < minDelta)
						delta += 2.0 * Math.PI;

					if (delta > maxDelta) {
						double distanceMax = delta - maxDelta;
						double distanceMin = minDelta - (delta - 2.0 * Math.PI);
						if (distanceMin < distanceMax) {
							delta = minDelta;
						} else {
							delta = maxDelta;
						}
					}
					rotateHelix(center, h.x, h.y, delta);

					// Re-assigns unpaired atoms
					boolean over = false;
					helixStart = _vp.getRNA().getCoords(h.x);
					helixStop = _vp.getRNA().getCoords(h.y);
					if (isDirect) {
						pHelR = computeAngle(center, helixStop) - base;
						pHelL = computeAngle(center, helixStart) - base;
					} else {
						pHelL = computeAngle(center, helixStop) - base;
						pHelR = computeAngle(center, helixStart) - base;
					}

					i = h.x - 1;
					Vector<Integer> nextBases = new Vector<Integer>();
					while (!over) {
						if (i < 0) {
							over = true;
						} else {
							if (_vp.getRNA().get_listeBases().get(i)
									.getElementStructure() == -1) {
								nextBases.add(new Integer(i));
							} else {
								over = true;
							}
						}
						i--;
					}
					Vector<Integer> prevBases = new Vector<Integer>();
					over = false;
					i = h.y + 1;
					while (!over) {
						if (i >= _vp.getRNA().get_listeBases().size()) {
							over = true;
						} else {
							if (_vp.getRNA().get_listeBases().get(i)
									.getElementStructure() == -1) {
								prevBases.add(new Integer(i));
							} else {
								over = true;
							}
						}
						i++;
					}
					double radius = center.distance(helixStart);
					double anglePrev, angleNext;
					double newAngle;
					if (isDirect) {
						anglePrev = normalizeAngle(pLimL - pHelR);
						angleNext = normalizeAngle(pHelL - pLimR);
					} else {
						anglePrev = normalizeAngle(pHelL - pLimR);
						angleNext = normalizeAngle(pLimL - pHelR);
					}

					for (i = 0; i < prevBases.size(); i++) {
						int k = prevBases.get(i).intValue();
						if (isDirect) {
							newAngle = base + pHelR + ((i + 1) * anglePrev)
									/ (prevBases.size() + 1);
						} else {
							newAngle = base + pHelL - ((i + 1) * anglePrev)
									/ (prevBases.size() + 1);
						}

						double newX = center.x + radius * Math.cos(newAngle);
						double newY = center.y + radius * Math.sin(newAngle);
						_vp.getRNA().setCoord(k, newX, newY);
					}
					for (i = 0; i < nextBases.size(); i++) {
						int k = nextBases.get(i).intValue();
						if (isDirect) {
							newAngle = base + pHelL - ((i + 1) * angleNext)
									/ (nextBases.size() + 1);
						} else {
							newAngle = base + pHelR + ((i + 1) * angleNext)
									/ (nextBases.size() + 1);
						}

						double newX = center.x + radius * Math.cos(newAngle);
						double newY = center.y + radius * Math.sin(newAngle);
						_vp.getRNA().setCoord(k, newX, newY);
					}
				}
			}
		}
	}


	private double normalizeAngle(double angle) {
		return normalizeAngle(angle, 0.0);
	}

	private double normalizeAngle(double angle, double base) {
		while (angle < base) {
			angle += 2.0 * Math.PI;
		}
		while (angle >= (2.0 * Math.PI) - base) {
			angle -= 2.0 * Math.PI;
		}
		return angle;
	}

	private double computeAngle(Point2D.Double center, Point2D.Double p) {
		double dist = center.distance(p);
		double angle = Math.asin((p.y - center.y) / dist);
		if (p.x - center.x < 0) {
			angle = Math.PI - angle;
		}
		return angle;
	}

	private void rotateHelix(Point2D.Double center, int i, int j, double angle) {
		for (int k = i; k <= j; k++) {
			Point2D.Double oldp = _vp.getRNA().getCoords(k);
			Point2D.Double newp = rotatePoint(center, oldp, angle);
			_vp.getRNA().setCoord(k, newp);
			if ((k != i) && (k != j)) {

				Point2D.Double oldc = _vp.getRNA().get_listeBases().get(k)
						.getCenter();
				Point2D.Double newc = rotatePoint(center, oldc, angle);
				_vp.getRNA().get_listeBases().get(k).setCenter(newc);
			}
		}
	}

	private Point2D.Double rotatePoint(Point2D.Double center, Point2D.Double p,
			double angle) {
		double dist = p.distance(center);
		double oldAngle = Math.asin((p.y - center.y) / dist);

		if (p.x - center.x < 0) {
			oldAngle = Math.PI - oldAngle;
		}

		double newX = (center.x + dist * Math.cos(oldAngle + angle));
		double newY = (center.y + dist * Math.sin(oldAngle + angle));

		return new Point2D.Double(newX, newY);
	}

	public void popupMenuCanceled(PopupMenuEvent arg0) {
	}

	public void popupMenuWillBecomeInvisible(PopupMenuEvent arg0) {
		//_vp.clearSelection();
		_vp.resetAnnotationHighlight();
		_vp.setSelectedBase(-1);
	}

	public void popupMenuWillBecomeVisible(PopupMenuEvent arg0) {
	}
}