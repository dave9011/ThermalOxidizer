
package main;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;

public abstract class AbstractValidator extends InputVerifier implements KeyListener {
    
    private JDialog popup;
    private Object parent;
    private JLabel messageLabel;
    private JLabel image;
    private Point point;
    private Dimension cDim;
    private Color color;
    private JTextField targetField = null;
    
    private boolean valid = true;
    
    private AbstractValidator(){
        color = new Color(243, 255, 159);
    }
    
    private AbstractValidator(JComponent c, String message) {
        this();  
        c.addKeyListener(AbstractValidator.this);
        messageLabel = new JLabel(message + " ");
        image = new JLabel(new ImageIcon(this.getClass().getResource("/img/error_icon.png")));
    }
    
    /**
     * @param parent A JDialog that implements the ValidationCapable interface.
     * @param c The JComponent to be validated.
     * @param message A message to be displayed in the pop up help tip if validation fails.
     */
    public AbstractValidator (JDialog parent, JComponent c, String message) {		
        this(c, message);
        this.parent = parent;
        popup = new JDialog(parent);
        initComponents();
    }
    
    /**
     * @param parent A JFrame that implements the ValidationCapable interface.
     * @param c The JComponent to be validated.
     * @param message A message to be displayed in the pop up help tip if validation fails.
     */
    public AbstractValidator (JFrame parent, JComponent c, String message) {
        this(c, message);
        this.parent = parent;
        popup = new JDialog(parent);
        initComponents();
    }
    
    protected abstract boolean validationCriteria(JComponent c);

    @Override
    public boolean verify(JComponent c) {
        
        if (!validationCriteria(c)) {
			
            if(parent instanceof WantsValidationStatus){ 
                ((WantsValidationStatus)parent).validateFailed();
            }
            
            c.setBackground(Color.PINK);
            c.setBorder(BorderFactory.createLineBorder(Color.RED));
            setPopup(c);
            popup.setVisible(true);
            return false;
            
        }
        
        if(targetField != null) {
            if(valid){
                    targetField.setBackground(Color.WHITE);  
                    targetField.setBorder( UIManager.getBorder("TextField.border") );
            }
            else{
                targetField.setBackground(Color.PINK); 
                targetField.setBorder(BorderFactory.createLineBorder(Color.RED));
            }
        }
    
        c.setBackground(Color.WHITE);
        c.setBorder( UIManager.getBorder("TextField.border") );
		
        if(parent instanceof WantsValidationStatus)
            ((WantsValidationStatus)parent).validatePassed();
		
        return true;
        
    }
    
    public void setValid(){   
        valid = true;
    }
    
    public void setNotValid(){    
        valid = false;
    }
    
    public void setTgtField(JTextField targetField){
        this.targetField = targetField;
    }
    
    private void setPopup(JComponent c){
        popup.setSize(0, 0);
        popup.setLocationRelativeTo(c);
        point = popup.getLocation();
        cDim = c.getSize();
        popup.setLocation(point.x-(int)cDim.getWidth()/2,
            point.y+(int)cDim.getHeight()/2);
        popup.pack();
    }
    
    /**
     * Changes the message that appears in the pop up help tip when a component's
     * data is invalid. Subclasses can use this to provide context sensitive help
     * depending on what the user did wrong.
     * 
     * @param message
     */
    protected void setMessage(String message) {
       messageLabel.setText(message);
    }

    @Override
    public void keyTyped(KeyEvent ke) {
        
    }

    @Override
    public void keyPressed(KeyEvent ke) {
        popup.setVisible(false);
    }

    @Override
    public void keyReleased(KeyEvent ke) {
       
    }
    
    private void initComponents() {
        popup.getContentPane().setLayout(new FlowLayout());
        popup.setUndecorated(true);
        popup.getContentPane().setBackground(color);
        popup.getContentPane().add(image);
        popup.getContentPane().add(messageLabel);
        popup.setFocusableWindowState(false);
    }
    
}
