package fr.orsay.lri.varna.applications.templateEditor;
import java.awt.BorderLayout;
import java.awt.ComponentOrientation;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.io.File;
import java.util.ArrayList;

import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.Icon;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JToggleButton;
import javax.swing.JToolBar;
import javax.swing.UIManager;
import javax.swing.filechooser.FileFilter;
import javax.swing.undo.UndoManager;

import fr.orsay.lri.varna.VARNAPanel;
import fr.orsay.lri.varna.applications.FileNameExtensionFilter;
import fr.orsay.lri.varna.exceptions.ExceptionInvalidRNATemplate;
import fr.orsay.lri.varna.exceptions.ExceptionNonEqualLength;
import fr.orsay.lri.varna.exceptions.ExceptionXMLGeneration;
import fr.orsay.lri.varna.models.templates.RNATemplateDrawingAlgorithmException;




public class TemplateEditor extends JFrame implements KeyListener, ActionListener {

	private TemplatePanel _sk;
	private VARNAPanel _vp;
	
	
	public TemplateEditor()
	{
		init();
	}
	
	private UndoManager manager;
	
	
	private void init()
	{
		try {
			_vp = new VARNAPanel();
		} catch (ExceptionNonEqualLength e) {
			e.printStackTrace();
		}
		JPanel p = new JPanel();
		p.setLayout(new GridLayout(1,2));
		
		JToolBar systemBar = new JToolBar();
		JButton loadButton = new JButton("Open...",UIManager.getIcon("FileView.directoryIcon")); 
		loadButton.setActionCommand("open");
		loadButton.addActionListener(this);
		loadButton.addKeyListener(this);
		JButton saveButton = new JButton("Save...",UIManager.getIcon("FileView.floppyDriveIcon")); 
		saveButton.setActionCommand("save");
		saveButton.addActionListener(this);
		saveButton.addKeyListener(this);
		JButton undoButton = new JButton("Undo"); 
		undoButton.setActionCommand("undo");
		undoButton.addActionListener(this);
		undoButton.addKeyListener(this);
		JButton redoButton = new JButton("Redo"); 
		redoButton.setActionCommand("redo");
		redoButton.addActionListener(this);
		redoButton.addKeyListener(this);
		JButton applyButton = new JButton("Apply"); 
		applyButton.setActionCommand("apply");
		applyButton.addActionListener(this);
		applyButton.addKeyListener(this);
		
		JButton flipButton = new JButton("Flip"); 
		flipButton.setActionCommand("flip");
		flipButton.addActionListener(this);
		flipButton.addKeyListener(this);

		
		systemBar.add(loadButton);
		systemBar.add(saveButton);
		systemBar.addSeparator();
		systemBar.add(undoButton);
		systemBar.add(redoButton);
		systemBar.addSeparator();
		systemBar.add(flipButton);
		systemBar.addSeparator();
		systemBar.add(applyButton);
		systemBar.addKeyListener(this);

		JToolBar toolBar = new JToolBar();
		ButtonGroup bg = new ButtonGroup(); 
		toolBar.setOrientation(JToolBar.VERTICAL);
		//toolBar.setLayout(new BoxLayout(toolBar, BoxLayout.PAGE_AXIS));
		JToggleButton selectButton = new JToggleButton("Select"); 
		selectButton.setActionCommand("select");
		selectButton.addActionListener(this);
		selectButton.addKeyListener(this);
		JToggleButton helixButton = new JToggleButton("Helix"); 
		helixButton.setActionCommand("helix");
		helixButton.addActionListener(this);
		helixButton.addKeyListener(this);
		JToggleButton unpairedButton = new JToggleButton("Unpaired"); 
		unpairedButton.setActionCommand("unpaired");
		unpairedButton.addActionListener(this);
		unpairedButton.addKeyListener(this);

		bg.add(selectButton);
		bg.add(helixButton);
		bg.add(unpairedButton);

		toolBar.add(selectButton);
		toolBar.add(helixButton);
		toolBar.add(unpairedButton);
		systemBar.addKeyListener(this);

		
		this.setLayout(new BorderLayout());
		_sk = new TemplatePanel();
		_sk.setPreferredSize(new Dimension(600,600));
		manager = new UndoManager();
		manager.setLimit(2000);
	   _sk.addUndoableEditListener(manager);
	   _sk.addKeyListener(this);
		
		JScrollPane jp = new JScrollPane(_sk,JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
		p.add(jp);
		p.add(_vp);
		getContentPane().add(systemBar,BorderLayout.PAGE_START);
		getContentPane().add(toolBar,BorderLayout.WEST);
		getContentPane().add(p,BorderLayout.CENTER);
		this.addKeyListener(this);
		
		_sk.requestFocusInWindow();
	}
	
	
	/**
	 * 
	 */
	private static final long serialVersionUID = -5942729520783690050L;

	public static void main(String[] argv)
	{
		TemplateEditor frame = new TemplateEditor();
		frame.pack();
		frame.setVisible(true);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	}

	public void keyPressed(KeyEvent e) {
		System.out.println(e);
		switch (e.getKeyCode())
		{
		  case (KeyEvent.VK_DELETE):
		  {
			  GraphicalTemplateElement h = _sk.getSelected();
			  _sk.Unselect();
			  _sk.getTemplateUI().removeElementUI(h);
			  _sk.repaint();
		  }
		  break;
		  case (KeyEvent.VK_Z):
		  {
			  if (e.isControlDown())
			  {
				  undo();
			  }
		  }
		  break;
		  case (KeyEvent.VK_Y):
		  {
			  if (e.isControlDown())
			  {
				  redo();
			  }
		  }
		  break;
		}
	}
	
	public void undo()
	{
		if (manager.canUndo())
		{
			System.out.println("Undo: "+manager.getUndoPresentationName());
			manager.undo();
		}
	}

	public void redo()
	{
		  if (manager.canRedo())
		  {
			  System.out.println("Redo: "+manager.getRedoPresentationName());
		      manager.redo();
		  }
	}	
	
	public void keyReleased(KeyEvent e) {
		System.out.println(e);
	}

	public void keyTyped(KeyEvent e) {
		System.out.println(e);
	}

	public void actionPerformed(ActionEvent e) {
		if (e.getActionCommand().equals("undo"))
		{
			undo();
		}
		else if (e.getActionCommand().equals("redo"))
		{
			redo();
		}
		else if (e.getActionCommand().equals("flip"))
		{
			GraphicalTemplateElement gr = _sk.getSelected();
			if (gr != null)
			{
				if (gr instanceof Helix)
				{
					_sk.getTemplateUI().flipHelixUI((Helix)gr);
				}
			}
		}
		else if (e.getActionCommand().equals("apply"))
		{
			try {
				_vp.getRNA().drawRNATemplate(_sk.getTemplate());
				_vp.repaint();
			} catch (RNATemplateDrawingAlgorithmException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}
		}
		else if (e.getActionCommand().equals("save"))
		{
			JFileChooser chooser = new JFileChooser();
			FileFilter filter = new FileNameExtensionFilter("VARNA RNA template (.xml)", "xml");
		    chooser.setFileFilter(filter);
			if (chooser.showSaveDialog(_sk) == JFileChooser.APPROVE_OPTION) {
				String path = chooser.getSelectedFile().getAbsolutePath();
				if (!path.toLowerCase().endsWith(".xml")) {
					path = path + ".xml";
				}
				try {
					_sk.getTemplate().toXMLFile(new File(path));
					System.out.println("Template saved in " + path);
				} catch (ExceptionXMLGeneration e1) {
					e1.printStackTrace();
				} catch (ExceptionInvalidRNATemplate e1) {
					e1.printStackTrace();
				}
			}
		}
		else if (e.getActionCommand().equals("open"))
		{
			JFileChooser chooser = new JFileChooser();
		    FileFilter filter = new FileNameExtensionFilter("VARNA RNA template (.xml)", "xml");
		    chooser.setFileFilter(filter);
			if (chooser.showOpenDialog(_sk) == JFileChooser.APPROVE_OPTION) {
				File templatePath = chooser.getSelectedFile();
				_sk.loadFromXmlFile(templatePath);
				
				// Empty the cancel history
				manager.discardAllEdits();
			}
		}
		
	}
}
