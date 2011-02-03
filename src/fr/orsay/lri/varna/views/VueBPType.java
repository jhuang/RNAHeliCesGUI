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
package fr.orsay.lri.varna.views;

import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;

import fr.orsay.lri.varna.VARNAPanel;
import fr.orsay.lri.varna.models.rna.ModeleStyleBP;


public class VueBPType implements ActionListener {

	private VARNAPanel _vp;
	ModeleStyleBP _msbp;
	private JComboBox _edge5;
	private JComboBox _edge3;
	private JComboBox _stericity;
	private JPanel panel;

	public VueBPType(VARNAPanel vp, ModeleStyleBP msbp) {
		_vp = vp;
		_msbp = msbp;

		ModeleStyleBP.Edge[] edges = ModeleStyleBP.Edge.values();
		ModeleStyleBP.Edge bck = msbp.getEdgePartner5();
		_edge5 = new JComboBox(edges);
		for (int i = 0; i < edges.length; i++) {
			if (edges[i] == bck)
				_edge5.setSelectedIndex(i);
		}

		bck = msbp.getEdgePartner3();
		_edge3 = new JComboBox(edges);
		for (int i = 0; i < edges.length; i++) {
			if (edges[i] == bck)
				_edge3.setSelectedIndex(i);
		}

		ModeleStyleBP.Stericity[] sters = ModeleStyleBP.Stericity.values();
		ModeleStyleBP.Stericity bcks = msbp.getStericity();
		_stericity = new JComboBox(sters);
		for (int i = 0; i < sters.length; i++) {
			if (sters[i] == bcks)
				_stericity.setSelectedIndex(i);
		}

		_edge5.addActionListener(this);
		_edge3.addActionListener(this);
		_stericity.addActionListener(this);

		panel = new JPanel();
		panel.setLayout(new FlowLayout(FlowLayout.LEFT));

		JLabel label5 = new JLabel("5' edge: ");
		JLabel label3 = new JLabel("3' edge: ");
		JLabel labelS = new JLabel("Stericity: ");

		panel.add(label5);
		panel.add(_edge5);
		panel.add(label3);
		panel.add(_edge3);
		panel.add(labelS);
		panel.add(_stericity);
	}

	public JPanel getPanel() {
		return panel;
	}

	public ModeleStyleBP.Edge getEdge5() {
		return (ModeleStyleBP.Edge) _edge5.getSelectedItem();
	}

	public ModeleStyleBP.Edge getEdge3() {
		return (ModeleStyleBP.Edge) _edge3.getSelectedItem();
	}

	public ModeleStyleBP.Stericity getStericity() {
		return (ModeleStyleBP.Stericity) _stericity.getSelectedItem();
	}

	public VARNAPanel get_vp() {
		return _vp;
	}

	public void actionPerformed(ActionEvent e) {
		_msbp.setEdge5(getEdge5());
		_msbp.setEdge3(getEdge3());
		_msbp.setStericity(getStericity());
		_vp.repaint();
	}

}
